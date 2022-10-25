#include <omp.h>
#include <spdlog/spdlog.h>

#include <Eigen/Dense>
#include <cxxopts.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "XYZGrid.h"
#include "gdal.h"
#include "mystructs.h"
#include "myutils.h"
#include "qin_io.h"

#define FIND_MEDIAN 0

namespace fs = std::filesystem;
namespace logg = spdlog;
using namespace std;
using namespace motionavg::TranslateND;

bool read_correspondence_pts_pix_binary_num(const char* filename, unsigned int& numpts) {
  std::fstream fs(filename, std::ios::in | std::ios::binary);
  if (!fs.is_open()) {
    std::cout << "file opening error" << std::endl;
    return false;
  }
  fs.read((char*)&numpts, sizeof(unsigned int));
  fs.close();
}
//
bool read_correspondence_pts_pix_binary(const char* filename, float* datapointerout) {
  std::fstream fs(filename, std::ios::in | std::ios::binary);
  if (!fs.is_open()) {
    std::cout << "file opening error" << std::endl;
    return false;
  }

  unsigned int numpt;
  fs.read((char*)&numpt, sizeof(unsigned int));
  fs.read((char*)datapointerout, sizeof(float) * numpt * 4);
  fs.close();

  return true;
}

bool write_correspondence_offset_binary(const char* filename, float* data, unsigned int num_pts) {
  std::ofstream ofs(filename, std::ios::binary);
  if (!ofs.is_open()) {
    std::cerr << "file open error " << filename << std::endl;
    return false;
  }

  ofs.write((char*)&num_pts, sizeof(decltype(num_pts)));
  ofs.write((char*)data, sizeof(float) * num_pts * 4);
  ofs.close();
  return true;
}

bool read_correspondence_offset_binary_num(const char* filename, unsigned int& numpts) {
  std::fstream fs(filename, std::ios::in | std::ios::binary);
  if (!fs.is_open()) {
    std::cout << "file opening error" << std::endl;
    return false;
  }
  fs.read((char*)&numpts, sizeof(unsigned int));
  fs.close();
}
//
bool read_correspondence_offset_binary(const char* filename, float* datapointerout) {
  std::fstream fs(filename, std::ios::in | std::ios::binary);
  if (!fs.is_open()) {
    std::cout << "file opening error" << std::endl;
    return false;
  }

  unsigned int numpt;
  fs.read((char*)&numpt, sizeof(unsigned int));
  fs.read((char*)datapointerout, sizeof(float) * numpt * 4);
  fs.close();

  return true;
}

vector<QinPose> readQinv2File(std::string input_qin) {
  ifstream ifs(input_qin);
  size_t num_pose;
  ifs >> num_pose;
  vector<QinPose> poses(num_pose);
  for (size_t i = 0; i < num_pose; ++i) ifs >> poses[i];
  return poses;
}

cxxopts::Options parseOptions(std::string exepath = "") {
  std::string exename = fs::path(exepath).filename().string();
  cxxopts::Options options(exename, "create_graph from folder");
  // clang-format off
  options.add_options()
    ("i,input_folder", "input folder", cxxopts::value<std::string>())
	("qinv2file", "qinv2 file", cxxopts::value<std::string>())
    ("gcpfile", "gcp json file", cxxopts::value<std::string>()->default_value(""))
	("o,output_name", "output filename", cxxopts::value<std::string>()->default_value("graph.txt"))
	("minrays", "minimum num of rays", cxxopts::value<int>()->default_value("3"))
	("minmatches", "minimum num of matches for valid pair", cxxopts::value<int>()->default_value("100"))
    ("sigma", "sigma of gaussian kernel for confidence [meter]", cxxopts::value<float>()->default_value("1"))
	("h,help", "Print help")
	("v,verbose", "verbose level (trace - 0, debug - 1, info - 2, warn - 3, error - 4, critical - 5, off - 6)", cxxopts::value<int>()->default_value("2"))
	;
  // clang-format on

  options.parse_positional({"input_folder", "qinv2file", "output_name"});
  options.positional_help("input_folder qinv2file output_name [options]");
  return options;
}

using GraphT = TranslateGraph3WithGCP;
int main(int argc, char** argv) {
  GDALAllRegister();

  float INVALID_VALUE = std::numeric_limits<float>::infinity();
  cxxopts::Options options = parseOptions(argv[0]);
  cxxopts::ParseResult args = options.parse(argc, argv);
  spdlog::set_level(static_cast<logg::level::level_enum>(args["verbose"].as<int>()));

  if (args.count("help") != 0) {
    cout << options.help() << endl;
    return 0;
  }
  if (args.count("input_folder") == 0 || args.count("qinv2file") == 0) {
    cout << options.help() << endl;
    return 0;
  }
  bool hasGCP = args.count("gcpfile") > 0;
  fs::path rootdirpath(args["input_folder"].as<string>());
  fs::path qinv2path(args["qinv2file"].as<string>());
  fs::path gcppath(args["gcpfile"].as<string>());
  fs::path outfilename(args["output_name"].as<string>()), outfilepath;
  int minimum_rays = args["minrays"].as<int>();
  bool preloadgrids = false;
  bool minimum_matches = args["minmatches"].as<int>();
  float sigma = args["sigma"].as<float>();
  float invsigmasq_2 = 2.f / std::pow(sigma, 2.f);
  if (outfilename.is_absolute())
    outfilepath = outfilename;
  else
    outfilepath = rootdirpath / outfilename;

  if (!fs::exists(qinv2path)) {
    spdlog::critical("QinPose Not found: {}", qinv2path.string());
    return 1;
  }

  if (hasGCP && !fs::exists(gcppath)) {
    spdlog::critical("GCP file Not Found: {}", gcppath.string());
    return 1;
  }

  fs::path lasdirpath = rootdirpath / "Point_clouds" / "point_clouds";
  fs::path tempfolderdirpath = rootdirpath / "Point_clouds" / "temp_folder_cluster";
  fs::path Neighborhood_file_path = lasdirpath / "Neighborhood_file.bin";
  if (!fs::exists(Neighborhood_file_path)) {
    spdlog::critical("File {} not found!", Neighborhood_file_path.string());
    return 1;
  }

  // Build graph
  // using GraphT = TranslateGraph<3>;
  GraphT g;
  g.basepath = rootdirpath.string();
  g.baseID = -1;

  vector<QinPose> poses = readQinv2File(qinv2path.string());
  // Add Nodes
  for (size_t i = 0; i < poses.size(); ++i) {
    GraphT::Node n;
    n.name = poses[i].imgname;
    fs::path tmpfolder = tempfolderdirpath / fs::path(n.name).stem();
    if (!fs::exists(tmpfolder)) {
      spdlog::error("Tmp folder {} not found!", tmpfolder.string());
      return 1;
    }
    n.path = fs::relative(tmpfolder, rootdirpath).string();
    n.xfm[0] = poses[i].x;
    n.xfm[1] = poses[i].y;
    n.xfm[2] = poses[i].z;
    g.insertNode(n);
  }
  
  if (hasGCP) {
    ifstream ifs;
    // Add GCP
    using json = nlohmann::json;
    ifs.open(gcppath.string());
    json gcpfilejson = json::parse(ifs);
    json gcpjson = gcpfilejson["GCP"];
    for (auto& [gcpname, gcpval] : gcpjson.items()) {
      GCPGraph::GCPPoint gcppt;
      gcppt.name = gcpname;
      gcppt.x = gcpval["x"].get<double>();
      gcppt.y = gcpval["y"].get<double>();
      gcppt.z = gcpval["z"].get<double>();

      gcppt.ex = gcpval.value("ex", 999.f);
      gcppt.ey = gcpval.value("ey", 999.f);
      gcppt.ez = gcpval.value("ez", 999.f);

      int gcpid = g.gcps.size();
      g.gcps.push_back(gcppt);

      for (auto& [key, view] : gcpval["views"].items()) {
        // find viewid
        std::string viewname = key;
        auto nodeit = std::find_if(g.nodes.begin(), g.nodes.end(), [&viewname](const GraphT::Node& n) { return viewname == n.name; });
        int viewid = std::distance(g.nodes.begin(), nodeit);
        GCPGraph::GCPLink gcplink;
        gcplink.gcpid = gcpid;
        gcplink.viewid = viewid;
        gcplink.u = view["u"].get<double>();
        gcplink.v = view["v"].get<double>();
        g.gcplinks.push_back(gcplink);
      }
    }

    // Read GCP constraints
    for (auto& l : g.gcplinks) {
      auto& view = g.nodes[l.viewid];
      auto& gcppt = g.gcps[l.gcpid];

      string viewName = fs::path(view.name).stem().string();

      // clang-format off
      fs::path viewgridpath[3] = {
        tempfolderdirpath / viewName / fmt::format("{}_Xgrid.tif", viewName),
        tempfolderdirpath / viewName / fmt::format("{}_Ygrid.tif", viewName),
        tempfolderdirpath / viewName / fmt::format("{}_Zgrid.tif", viewName)};
      // clang-format on
      XYZGrid viewXYZ(false);
      viewXYZ.open(viewgridpath, fs::path());
      
      uint16_t dummy;
      viewXYZ.sample(l.u, l.v, l.dx, l.dy, l.dz, dummy);
      // TODO: The erro could be derived from grids and correspondences.
      l.ex = 0.01;
      l.ey = 0.01;
      l.ez = 0.01;
      double viewX = l.dx + view.xfm[0];
      double viewY = l.dy + view.xfm[1];
      double viewZ = l.dz + view.xfm[2];

      spdlog::trace("GCP[{}]-View[{}]: {} {} {}", gcppt.name, view.name, viewX - gcppt.x, viewY - gcppt.y, viewZ - gcppt.z);
    }
  }
  
  // Add Edges
  int numimg, maxpair;
  std::vector<int> neighborbin;
  {
    ifstream ifs;
    ifs.open(Neighborhood_file_path.string(), ios::binary);
    ifs.read(reinterpret_cast<char*>(&numimg), sizeof(int));
    ifs.read(reinterpret_cast<char*>(&maxpair), sizeof(int));
    spdlog::debug("neighbor bin: {} {}", numimg, maxpair);
    neighborbin.resize(numimg * (maxpair + 1));
    ifs.read(reinterpret_cast<char*>(neighborbin.data()), numimg * (maxpair + 1) * sizeof(int));
    ifs.close();
  }

  for (int i = 0; i < numimg; ++i)  // i-th image
  {
    const int* ids = &neighborbin[i * (maxpair + 1)];
    spdlog::trace("{}: {}", i, fmt::join(ids + 1, ids + maxpair + 1, ", "));

    int srcId = i;
    string srcName = fs::path(g.nodes[srcId].name).stem().string();

    // if (srcName != "029_065_id3188c81615_124544_Backward") continue;

    // clang-format off
    fs::path srcgridpath[3] = {
        tempfolderdirpath / srcName / fmt::format("{}_Xgrid.tif", srcName),
        tempfolderdirpath / srcName / fmt::format("{}_Ygrid.tif", srcName),
        tempfolderdirpath / srcName / fmt::format("{}_Zgrid.tif", srcName)};
    // clang-format on
    fs::path srcraynumpath = tempfolderdirpath / srcName / fmt::format("{}_ray_num.tif", srcName);

    for (int _d = 0; _d < 3; ++_d)
      if (!fs::exists(srcgridpath[_d])) {
        spdlog::error("File Not Found: {}", srcgridpath[_d].string());
        return 1;
      }
    if (!fs::exists(srcraynumpath)) {
      spdlog::error("File Not Found: {}", srcraynumpath.string());
      return 1;
    }

    XYZGrid srcXYZ(true);
    srcXYZ.open(srcgridpath, srcraynumpath);

#pragma omp parallel for
    for (int j = 1; j < maxpair + 1; ++j)  // j-th neighbor
    {
      if (ids[j] == -1) continue;
      int tgtId = ids[j];
      GraphT::Edge e;
      e.source = srcId;
      e.target = tgtId;
      string tgtName = fs::path(g.nodes[tgtId].name).stem().string();
      // if (tgtName != "029_056_id3179c82579_124523_Nadir") continue;
      fs::path correspondencepixbin = tempfolderdirpath / srcName / fmt::format("{}_{}_correspondence_pix.bin", srcName, tgtName);
      fs::path correspondenceoffsetbin = tempfolderdirpath / srcName / fmt::format("{}_{}_correspondence_offset.bin", srcName, tgtName);

      e.name = fs::relative(correspondencepixbin, rootdirpath).string();
      if (!fs::exists(correspondencepixbin)) {
        spdlog::error("Cannot find correspondence file: {}", correspondencepixbin.string());
        continue;
      }

      unsigned int numpts;
      read_correspondence_pts_pix_binary_num(correspondencepixbin.string().c_str(), numpts);
      spdlog::trace("numpts {}", numpts);
      if (numpts < minimum_matches) continue;
      std::unique_ptr<float> matches(new float[4 * numpts]);
      std::unique_ptr<float> offset_n_weight(new float[4 * numpts]);
      read_correspondence_pts_pix_binary(correspondencepixbin.string().c_str(), matches.get());

      spdlog::trace("src {} tgt {}", srcName, tgtName);

      // clang-format off
      fs::path tgtgridpath[3] = {
          tempfolderdirpath / tgtName / fmt::format("{}_Xgrid.tif", tgtName),
          tempfolderdirpath / tgtName / fmt::format("{}_Ygrid.tif", tgtName),
          tempfolderdirpath / tgtName / fmt::format("{}_Zgrid.tif", tgtName)};
      // clang-format on
      fs::path tgtraynumpath = tempfolderdirpath / tgtName / fmt::format("{}_ray_num.tif", tgtName);

      // clang-format off
      fs::path srctempgridpath[3] = {
          tempfolderdirpath / srcName / fmt::format("{}_{}_temgridX.tif", srcName, tgtName),
          tempfolderdirpath / srcName / fmt::format("{}_{}_temgridY.tif", srcName, tgtName),
          tempfolderdirpath / srcName / fmt::format("{}_{}_temgridZ.tif", srcName, tgtName)};
      // clang-format on

      for (int _d = 0; _d < 3; ++_d)
        if (!fs::exists(tgtgridpath[_d])) {
          spdlog::error("File Not Found: {}", tgtgridpath[_d].string());
          continue;
        }
      if (!fs::exists(tgtraynumpath)) {
        spdlog::error("File Not Found: {}", tgtraynumpath.string());
        continue;
      }

      for (int _d = 0; _d < 3; ++_d)
        if (!fs::exists(srctempgridpath[_d])) {
          spdlog::error("File Not Found: {}", srctempgridpath[_d].string());
          continue;
        }

      XYZGrid srcTmpXYZ(preloadgrids);
      srcTmpXYZ.open(srctempgridpath, fs::path());
      XYZGrid tgtXYZ(preloadgrids);
      tgtXYZ.open(tgtgridpath, tgtraynumpath);

      double meanSrc[3] = {0}, meanTgt[3] = {0};
#if FIND_MEDIAN
      std::vector<float> offset[3], offDist;
      offset[0].reserve(numpts);
      offset[1].reserve(numpts);
      offset[2].reserve(numpts);
      offDist.reserve(numpts);
#endif

      int numMatch = 0;
      float totalWeight = 0;
      for (int k = 0; k < numpts; k += 1)  // k-th matched points
      {
        std::fill_n(offset_n_weight.get() + 4 * k, 3, INVALID_VALUE);
        offset_n_weight.get()[4 * k + 3] = 0;
        float& srcX = matches.get()[4 * k + 0];
        float& srcY = matches.get()[4 * k + 1];
        float& tgtX = matches.get()[4 * k + 2];
        float& tgtY = matches.get()[4 * k + 3];

        float srcMatch[3], srcTmpMatch[3], tgtMatch[3];
        uint16_t srcNumRay, srcTmpNumRay, tgtNumRay;
        bool valid0 = srcXYZ.sample(srcX, srcY, srcMatch[0], srcMatch[1], srcMatch[2], srcNumRay);
        if (!valid0) continue;
        bool valid1 = tgtXYZ.sample(tgtX, tgtY, tgtMatch[0], tgtMatch[1], tgtMatch[2], tgtNumRay);
        if (!valid1) continue;
        bool valid00 = srcTmpXYZ.sample(srcX, srcY, srcTmpMatch[0], srcTmpMatch[1], srcTmpMatch[2], srcTmpNumRay);
        if (!valid00) continue;

        if (srcNumRay < minimum_rays || tgtNumRay < minimum_rays)
          continue;
        else {
          ++numMatch;
          float _weight = std::pow(srcTmpMatch[0] - srcMatch[0], 2) + std::pow(srcTmpMatch[1] - srcMatch[1], 2) + std::pow(srcTmpMatch[2] - srcMatch[2], 2);
          _weight = std::exp(-invsigmasq_2 * _weight);
          totalWeight += _weight;

          offset_n_weight.get()[4 * k + 0] = srcMatch[0] - tgtMatch[0];
          offset_n_weight.get()[4 * k + 1] = srcMatch[1] - tgtMatch[1];
          offset_n_weight.get()[4 * k + 2] = srcMatch[2] - tgtMatch[2];
          offset_n_weight.get()[4 * k + 3] = _weight;
#if FIND_MEDIAN
          double _distsq = 0;
#endif
          for (int _d = 0; _d < 3; ++_d) {
            meanSrc[_d] += _weight * srcMatch[_d];
            meanTgt[_d] += _weight * tgtMatch[_d];
#if FIND_MEDIAN
            offset[_d].emplace_back(srcMatch[_d] - tgtMatch[_d]);
            _distsq += std::pow(srcMatch[_d] - tgtMatch[_d], 2);
          }
          offDist.emplace_back(_distsq);
#else
          }
#endif
        }
      }

      if (numMatch < minimum_matches) continue;
      bool valid = true;

#if FIND_MEDIAN
      // find median value
      std::vector<float> offDistCopy = offDist;
      auto middle = offDistCopy.begin() + (offDistCopy.size() / 2);
      std::nth_element(offDistCopy.begin(), middle, offDistCopy.end());
      float nthvalue = *middle;
      auto it = std::find(offDist.begin(), offDist.end(), nthvalue);
      auto pos = std::distance(offDist.begin(), it);

      double medianXfm[3] = {offset[0][pos], offset[1][pos], offset[2][pos]};

      spdlog::trace("Median: {} {} {}", medianXfm[0], medianXfm[1], medianXfm[2]);
#endif
      // compute mean value
      for (int _d = 0; _d < 3; ++_d) {
        meanSrc[_d] /= totalWeight;
        meanTgt[_d] /= totalWeight;
        e.xfm[_d] = meanSrc[_d] - meanTgt[_d];
        e.cov[4 * _d] = 1. / totalWeight;
        if (std::isnan(e.xfm[_d]) || std::isinf(e.xfm[_d])) valid = false;
      }
#if FIND_MEDIAN
      spdlog::trace("Mean: {} {} {}", meanSrc[0] - meanTgt[0], meanSrc[1] - meanTgt[1], meanSrc[2] - meanTgt[2]);
#endif
      if (!valid) continue;

      double initXfm[3] = {g.nodes[tgtId].xfm[0] - g.nodes[srcId].xfm[0], g.nodes[tgtId].xfm[1] - g.nodes[srcId].xfm[1], g.nodes[tgtId].xfm[2] - g.nodes[srcId].xfm[2]};

      double diffXfm[3] = {std::abs(initXfm[0] - e.xfm[0]), std::abs(initXfm[1] - e.xfm[1]), std::abs(initXfm[2] - e.xfm[2])};
      spdlog::trace("Init: {} {} {}", initXfm[0], initXfm[1], initXfm[2]);
      spdlog::trace("Edge: {} {} {}", e.xfm[0], e.xfm[1], e.xfm[2]);
      spdlog::debug("Diff: {} {} {} cm {}, {}, {}", diffXfm[0] * 100, diffXfm[1] * 100, diffXfm[2] * 100, srcName, tgtName, totalWeight);

#pragma omp critical
      g.insertEdge(e);

      write_correspondence_offset_binary(correspondenceoffsetbin.string().c_str(), offset_n_weight.get(), numpts);
    }
  }

  ofstream ofs(outfilepath.string());
  ofs << &g;
  ofs.close();
  spdlog::info("Done");

  return 0;
}