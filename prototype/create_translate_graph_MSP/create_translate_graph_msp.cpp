// #include <omp.h>
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

template< typename T>
class StatVar {

std::string _name;
  T _E;
  T _sqE;
  uint32_t _cnt;
 public:
  StatVar(std::string name):_name(name) { _cnt = 0; }
  void measure(T _v) {
    if (_cnt == 0) {
      _E = _v;
      _sqE = _v.cwiseProduct(_v);
    } else {
      _E += (_v - _E) / (_cnt + 1);
      _sqE += (_v.cwiseProduct(_v) - _sqE) / (_cnt + 1);
    }
    _cnt += 1;
  }

  T getE() const { return _E; }
  T getVar() const { return _sqE - _E.cwiseProduct(_E); }
  uint32_t getCount() const { return _cnt; }
  
  void print(std::ostream& out) {
    out << "E(" << _name << ")=" << getE().transpose() <<
        ", D(" << _name << ")=" << getVar().transpose() << " of " << getCount() << " measurements" << std::endl;
  }

  void print(spdlog::logger& logger) { 
    std::stringstream ss;
    ss << "E(" << _name << ")=" << getE().transpose() << ", D(" << _name << ")=" << getVar().transpose() << " of " << getCount() << " measurements";
    logger.info(ss.str());
  }
};

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
  fs.read((char*)datapointerout, sizeof(float) * numpt * 7);
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

float INVALID_VALUE = std::numeric_limits<float>::infinity();
/*
bool create_edge_minimize_reproject_error(const XYZGrid* srcXYZ, const XYZGrid* srcTmpXYZ, const XYZGrid* tgtTmpXYZ, const XYZGrid* tgtXYZ, const QinPose& srcPose, const QinPose& tgtPose,
                                          const float* matches, int numpts, float* offset_n_weight, int minimum_rays, float invsigmasq_2, int minimum_matches, double* xfm, double* cov) {
  auto srcK = srcPose.GetK();
  auto srcR = srcPose.GetR();
  auto srcC = srcPose.GetC();
  auto tgtK = tgtPose.GetK();
  auto tgtR = tgtPose.GetR();
  auto tgtC = tgtPose.GetC();

  Eigen::Vector3d init_bij = -srcC + tgtC;
  Eigen::Vector3d expected_bij{0, 0, 0};
  int num_measure = 0;
  float totalWeight = 0;
  for (int k = 0; k < numpts; k += 1)  // k-th matched points
  {
    std::fill_n(offset_n_weight + 4 * k, 3, INVALID_VALUE);
    offset_n_weight[4 * k + 3] = 0;
    const float& srcX = matches[7 * k + 0];
    const float& srcY = matches[7 * k + 1];
    const float& tgtX = matches[7 * k + 2];
    const float& tgtY = matches[7 * k + 3];

    float srcMatch[3];
    // float srcTmpMatch[3];
    // float tgtMatch[3];
    // float tgtTmpMatch[3];
    uint16_t srcNumRay;
    // uint16_t srcTmpNumRay;
    // uint16_t tgtNumRay;
    // uint16_t tgtTmpNumRay;
    bool valid0 = srcXYZ->sample(srcX, srcY, srcMatch[0], srcMatch[1], srcMatch[2], srcNumRay);
    if (!valid0) continue;
    //bool valid1 = tgtXYZ->sample(tgtX, tgtY, tgtMatch[0], tgtMatch[1], tgtMatch[2], tgtNumRay);
    //if (!valid1) continue;
    //bool valid00 = srcTmpXYZ->sample(srcX, srcY, srcTmpMatch[0], srcTmpMatch[1], srcTmpMatch[2], srcTmpNumRay);
    //if (!valid00) continue;

    if (srcNumRay < minimum_rays
        //|| tgtNumRay < minimum_rays
    )
      continue;
    else {
      Eigen::Vector3d init_bij = -srcC + tgtC;
      Eigen::Vector3d objpt_src = Eigen::Vector3d(srcMatch[0], srcMatch[1], srcMatch[2]);
      // Eigen::Vector3d objpt_tgt = Eigen::Vector3d(tgtMatch[0], tgtMatch[1], tgtMatch[2]);
      // Eigen::Vector3d campt_tgt2tgt = tgtR * objpt_tgt;
      // Eigen::Vector3d campt_tmptgt = tgtR * Eigen::Vector3d(tgtTmpMatch[0], tgtTmpMatch[1], tgtTmpMatch[2]);
      Eigen::Vector3d campt_src2tgt = tgtR * (srcC + objpt_src - tgtC);
      // Eigen::Vector3d campt_srctmp2tgt = tgtR * (srcC + Eigen::Vector3d(srcTmpMatch[0], srcTmpMatch[1], srcTmpMatch[2]) - tgtC);
      double depth_tgt2tgt = campt_tgt2tgt[2];
      double depth_src2tgt = campt_src2tgt[2];
      // Eigen::Vector3d imgpt_tgt2tgt = tgtK * campt_tgt2tgt / depth_tgt2tgt;
      Eigen::Vector3d imgpt_src2tgt = tgtK * campt_src2tgt / depth_src2tgt;

      Eigen::Vector3d _bij = objpt_src - tgtR.transpose() * tgtK.inverse() * Eigen::Vector3d{depth_src2tgt * tgtX, depth_src2tgt * tgtY, depth_src2tgt};

      if (num_measure == 0) {
        expected_bij = _bij;
      } else {
        expected_bij += (_bij - expected_bij) / (num_measure + 1);
      }

      ++num_measure;
    }
  }

  if (num_measure < minimum_matches) return false;
  bool valid = true;

  totalWeight = num_measure;
  for (int _d = 0; _d < 3; ++_d) {
    xfm[_d] = expected_bij[_d];
    cov[4 * _d] = 1. / totalWeight;
    if (std::isnan(xfm[_d]) || std::isinf(xfm[_d])) valid = false;
  }

  return valid;
}*/

StatVar<Eigen::Vector2d> total_stat_reprojection_error("total reproj");

bool create_edge_minimize_reproject_error2(const XYZGrid* srcTmpXYZ, const XYZGrid* tgtTmpXYZ, const QinPose& srcPose, const QinPose& tgtPose, const float* matches, int numpts, float* offset_n_weight,
                                           float invsigmasq_2, int minimum_matches, double* xfm, double* cov, spdlog::logger* logger=nullptr) {
  auto srcK = srcPose.GetK();
  auto srcR = srcPose.GetR();
  auto srcC = srcPose.GetC();
  auto tgtK = tgtPose.GetK();
  auto tgtR = tgtPose.GetR();
  auto tgtC = tgtPose.GetC();

  Eigen::Vector3d init_bij = -srcC + tgtC;
  Eigen::Vector3d expected_bij{0, 0, 0};
  int num_measure = 0;
  float totalWeight = 0;
  
  StatVar<Eigen::Vector2d> stat_src2tgt_error("src2tgt_reproj");
  StatVar<Eigen::Vector2d> stat_tgt2tgt_error("tgt2tgt_reproj");
  StatVar<Eigen::Vector2d> stat_cortgt_error("cortgt_reproj");
  StatVar<Eigen::Vector3d> stat_bij("bij");

  for (int k = 0; k < numpts; k ++)  // k-th matched points
  {
    std::fill_n(offset_n_weight + 4 * k, 3, INVALID_VALUE);
    offset_n_weight[4 * k + 3] = 0;
    const float& srcX = matches[7 * k + 0];
    const float& srcY = matches[7 * k + 1];
    const float& tgtX = matches[7 * k + 2];
    const float& tgtY = matches[7 * k + 3];

    float srcTmpMatch[3];
    float tgtTmpMatch[3];
    uint16_t dummy;
    bool valid0 = srcTmpXYZ->sample(srcX, srcY, srcTmpMatch[0], srcTmpMatch[1], srcTmpMatch[2], dummy);
    if (!valid0) continue;
    bool valid1 = tgtTmpXYZ->sample(tgtX, tgtY, tgtTmpMatch[0], tgtTmpMatch[1], tgtTmpMatch[2], dummy);
    if (!valid1) continue;

    Eigen::Vector3d init_bij = -srcC + tgtC;
    Eigen::Vector3d objpt_src = Eigen::Vector3d(srcTmpMatch[0], srcTmpMatch[1], srcTmpMatch[2]);
    Eigen::Vector3d objpt_tgt = Eigen::Vector3d(tgtTmpMatch[0], tgtTmpMatch[1], tgtTmpMatch[2]);
    Eigen::Vector3d campt_src2tgt = tgtR * (srcC + objpt_src - tgtC);
    Eigen::Vector3d campt_tgt2tgt = tgtR * objpt_tgt;

    double depth_src2tgt = campt_src2tgt[2];
    double depth_tgt2tgt = campt_tgt2tgt[2];

    Eigen::Vector3d imgpt_tgt2tgt = tgtK * campt_tgt2tgt / depth_tgt2tgt;
    Eigen::Vector3d imgpt_src2tgt = tgtK * campt_src2tgt / depth_src2tgt;

    

    double depth_ = depth_tgt2tgt;

    Eigen::Vector3d _bij = objpt_src - tgtR.transpose() * tgtK.inverse() * (depth_ * Eigen::Vector3d{tgtX, tgtY, 1.0});

    // trust 
    //std::cout << "A:"<< imgpt_tgt2tgt[0] - tgtX << "," << imgpt_tgt2tgt[1] - tgtY << std::endl;
    //std::cout << "B:" << imgpt_src2tgt[0] - tgtX << "," << imgpt_src2tgt[1] - tgtY << std::endl;
    
    Eigen::Vector3d imgpt_correct_tgt = tgtK * tgtR * (objpt_src - _bij);
    imgpt_correct_tgt /= imgpt_correct_tgt[2];
    //std::cout << "C:" << imgpt_correct_tgt[0] - tgtX << "," << imgpt_correct_tgt[1] - tgtY << std::endl;

    stat_src2tgt_error.measure(Eigen::Vector2d{imgpt_src2tgt[0] - tgtX, imgpt_src2tgt[1] - tgtY});
    total_stat_reprojection_error.measure(Eigen::Vector2d{imgpt_src2tgt[0] - tgtX, imgpt_src2tgt[1] - tgtY});
    stat_tgt2tgt_error.measure(Eigen::Vector2d{imgpt_tgt2tgt[0] - tgtX, imgpt_tgt2tgt[1] - tgtY});
    stat_cortgt_error.measure(Eigen::Vector2d{imgpt_correct_tgt[0] - tgtX, imgpt_correct_tgt[1] - tgtY});
    stat_bij.measure(_bij);
    if (num_measure == 0) {
      expected_bij = _bij;
    } else {
      expected_bij += (_bij - expected_bij) / (num_measure + 1);
    }
    ++num_measure;
  }
  
  if (logger == nullptr) {
    stat_bij.print(std::cout);
    stat_src2tgt_error.print(std::cout);
    stat_tgt2tgt_error.print(std::cout);
    stat_cortgt_error.print(std::cout);
  } else {
    stat_bij.print(*logger);
    stat_src2tgt_error.print(*logger);
    stat_tgt2tgt_error.print(*logger);
    stat_cortgt_error.print(*logger);
  }

  if (num_measure < minimum_matches) return false;
  bool valid = true;

  totalWeight = num_measure;
  std::fill_n(cov, 9, 0.);
  for (int _d = 0; _d < 3; ++_d) {
    xfm[_d] = expected_bij[_d];
    cov[4 * _d] = 1. / totalWeight;
    if (std::isnan(xfm[_d]) || std::isinf(xfm[_d])) valid = false;
  }

  return valid;
}

bool create_edge_minimize_object_points(const XYZGrid* srcXYZ, const XYZGrid* srcTmpXYZ, const XYZGrid* tgtTmpXYZ, const XYZGrid* tgtXYZ, const float* matches, int numpts, float* offset_n_weight,
                                        int minimum_rays, float invsigmasq_2, int minimum_matches, double* xfm, double* cov) {
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
    std::fill_n(offset_n_weight + 4 * k, 3, INVALID_VALUE);
    offset_n_weight[4 * k + 3] = 0;
    const float& srcX = matches[4 * k + 0];
    const float& srcY = matches[4 * k + 1];
    const float& tgtX = matches[4 * k + 2];
    const float& tgtY = matches[4 * k + 3];

    float srcMatch[3], srcTmpMatch[3], tgtMatch[3];
    uint16_t srcNumRay, srcTmpNumRay, tgtNumRay;
    bool valid0 = srcXYZ->sample(srcX, srcY, srcMatch[0], srcMatch[1], srcMatch[2], srcNumRay);
    if (!valid0) continue;
    bool valid1 = tgtXYZ->sample(tgtX, tgtY, tgtMatch[0], tgtMatch[1], tgtMatch[2], tgtNumRay);
    if (!valid1) continue;
    bool valid00 = srcTmpXYZ->sample(srcX, srcY, srcTmpMatch[0], srcTmpMatch[1], srcTmpMatch[2], srcTmpNumRay);
    if (!valid00) continue;

    if (srcNumRay < minimum_rays || tgtNumRay < minimum_rays)
      continue;
    else {
      ++numMatch;
      float _weight = std::pow(srcTmpMatch[0] - srcMatch[0], 2) + std::pow(srcTmpMatch[1] - srcMatch[1], 2) + std::pow(srcTmpMatch[2] - srcMatch[2], 2);
      _weight = std::exp(-invsigmasq_2 * _weight);
      totalWeight += _weight;

      offset_n_weight[4 * k + 0] = srcMatch[0] - tgtMatch[0];
      offset_n_weight[4 * k + 1] = srcMatch[1] - tgtMatch[1];
      offset_n_weight[4 * k + 2] = srcMatch[2] - tgtMatch[2];
      offset_n_weight[4 * k + 3] = _weight;
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

  if (numMatch < minimum_matches) return false;
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

  logger.trace("Median: {} {} {}", medianXfm[0], medianXfm[1], medianXfm[2]);
#endif
  // compute mean value
  for (int _d = 0; _d < 3; ++_d) {
    meanSrc[_d] /= totalWeight;
    meanTgt[_d] /= totalWeight;
    xfm[_d] = meanSrc[_d] - meanTgt[_d];
    cov[4 * _d] = 1. / totalWeight;
    if (std::isnan(xfm[_d]) || std::isinf(xfm[_d])) valid = false;
  }
#if FIND_MEDIAN
  logger.trace("Mean: {} {} {}", meanSrc[0] - meanTgt[0], meanSrc[1] - meanTgt[1], meanSrc[2] - meanTgt[2]);
#endif
  if (!valid) return false;

  return true;
}

using GraphT = TranslateGraph3WithGCP;

#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>


int main(int argc, char** argv) {
  GDALAllRegister();

  cxxopts::Options options = parseOptions(argv[0]);
  cxxopts::ParseResult args = options.parse(argc, argv);

  

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

  
  auto console_sink = std::make_shared<spdlog::sinks::stderr_color_sink_mt>();
  console_sink->set_level(spdlog::level::debug);

  auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(outfilepath.string() + ".create_graph.log");
  file_sink->set_level(spdlog::level::trace);

  spdlog::logger logger("LOG", {console_sink, file_sink});
  logger.set_level(static_cast<logg::level::level_enum>(args["verbose"].as<int>()));


  if (!fs::exists(qinv2path)) {
    logger.critical("QinPose Not found: {}", qinv2path.string());
    return 1;
  }

  if (hasGCP && !fs::exists(gcppath)) {
    logger.critical("GCP file Not Found: {}", gcppath.string());
    return 1;
  }

  fs::path lasdirpath = rootdirpath / "Point_clouds" / "point_clouds";
  fs::path tempfolderdirpath = rootdirpath / "Point_clouds" / "temp_folder_cluster";
  fs::path Neighborhood_file_path = lasdirpath / "Neighborhood_file.bin";
  if (!fs::exists(Neighborhood_file_path)) {
    logger.critical("File {} not found!", Neighborhood_file_path.string());
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
      logger.error("Tmp folder {} not found!", tmpfolder.string());
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

      logger.trace("GCP[{}]-View[{}]: {} {} {}", gcppt.name, view.name, viewX - gcppt.x, viewY - gcppt.y, viewZ - gcppt.z);
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
    logger.debug("neighbor bin: {} {}", numimg, maxpair);
    neighborbin.resize(numimg * (maxpair + 1));
    ifs.read(reinterpret_cast<char*>(neighborbin.data()), numimg * (maxpair + 1) * sizeof(int));
    ifs.close();
  }

  std::vector<std::vector<int>> decoded_neighbors(numimg);
  for (int i = 0; i < numimg; ++i)  // i-th image
  {
    const int* ids = &neighborbin[i * (maxpair + 1)];
    for (const int* _pj = ids + 1; _pj < ids + maxpair + 1; ++_pj) {
      if (*_pj == -1) break;
      decoded_neighbors[i].push_back(*_pj);
    }
  }

  for (int i = 0; i < numimg; ++i)  // i-th image
  {
    const int* ids = &neighborbin[i * (maxpair + 1)];
    int srcId = i;
    string srcName = fs::path(g.nodes[srcId].name).stem().string();

    for (int j : decoded_neighbors[i]) {
      int tgtId = j;
      string tgtName = fs::path(g.nodes[tgtId].name).stem().string();
      GraphT::Edge e;
      e.source = srcId;
      e.target = tgtId;

      fs::path correspondencepixbin = tempfolderdirpath / srcName / fmt::format("{}_{}_correspondence_pix.bin", srcName, tgtName);
      fs::path correspondenceoffsetbin = tempfolderdirpath / srcName / fmt::format("{}_{}_correspondence_offset.bin", srcName, tgtName);

      e.name = fs::relative(correspondencepixbin, rootdirpath).string();
      if (!fs::exists(correspondencepixbin)) {
        logger.error("Cannot find correspondence file: {}", correspondencepixbin.string());
        continue;
      }

      unsigned int numpts;
      read_correspondence_pts_pix_binary_num(correspondencepixbin.string().c_str(), numpts);
      logger.trace("numpts {}", numpts);
      if (numpts < minimum_matches) continue;
      std::unique_ptr<float> matches(new float[7 * numpts]);
      std::unique_ptr<float> offset_n_weight(new float[7 * numpts]);
      read_correspondence_pts_pix_binary(correspondencepixbin.string().c_str(), matches.get());

      // Now, make combination of src-neighbors and tgt-neighbors, except src-tgt
      std::vector<std::pair<int, int>> xyzgrid_pairs;
      for (int adj_srcId : decoded_neighbors[srcId])
        for (int adj_tgtId : decoded_neighbors[tgtId]) {
          if (adj_srcId != tgtId && adj_tgtId != srcId) xyzgrid_pairs.emplace_back(adj_srcId, adj_tgtId);
        }
      logger.info("Edge[{},{}] requires {} grid pairs", e.source, e.target, xyzgrid_pairs.size());

      Eigen::Vector3d all_xfm = Eigen::Vector3d::Zero();
      Eigen::Matrix3d all_cov = Eigen::Matrix3d::Identity();
      all_cov.fill(10000);
      for (auto& p : xyzgrid_pairs) {
        int adj_srcId, adj_tgtId;
        std::tie(adj_srcId, adj_tgtId) = p;
        std::string adj_srcName = fs::path(g.nodes[adj_srcId].name).stem().string();
        std::string adj_tgtName = fs::path(g.nodes[adj_tgtId].name).stem().string();
        // clang-format off
    #if 1
        fs::path srctempgridpath[3] = {
          tempfolderdirpath / srcName / fmt::format("{}_{}_temgridX.tif", srcName, adj_srcName),
          tempfolderdirpath / srcName / fmt::format("{}_{}_temgridY.tif", srcName, adj_srcName),
          tempfolderdirpath / srcName / fmt::format("{}_{}_temgridZ.tif", srcName, adj_srcName)};
        fs::path tgttempgridpath[3] = {
          tempfolderdirpath / tgtName / fmt::format("{}_{}_temgridX.tif", tgtName, adj_tgtName),
          tempfolderdirpath / tgtName / fmt::format("{}_{}_temgridY.tif", tgtName, adj_tgtName),
          tempfolderdirpath / tgtName / fmt::format("{}_{}_temgridZ.tif", tgtName, adj_tgtName)};
    #else
        fs::path srctempgridpath[3] = {
          tempfolderdirpath / srcName / fmt::format("{}_{}_temgridX.tif", srcName, tgtName),
          tempfolderdirpath / srcName / fmt::format("{}_{}_temgridY.tif", srcName, tgtName),
          tempfolderdirpath / srcName / fmt::format("{}_{}_temgridZ.tif", srcName, tgtName)};
        fs::path tgttempgridpath[3] = {
          tempfolderdirpath / tgtName / fmt::format("{}_{}_temgridX.tif", tgtName, srcName),
          tempfolderdirpath / tgtName / fmt::format("{}_{}_temgridY.tif", tgtName, srcName),
          tempfolderdirpath / tgtName / fmt::format("{}_{}_temgridZ.tif", tgtName, srcName)};
    #endif
        // clang-format on
        for (int _d = 0; _d < 3; ++_d) {
          if (!fs::exists(srctempgridpath[_d])) {
            logger.error("File Not Found: {}", srctempgridpath[_d].string());
            return 1;
          }
          if (!fs::exists(tgttempgridpath[_d])) {
            logger.error("File Not Found: {}", tgttempgridpath[_d].string());
            return 1;
          }
        }  // for (int _d = 0; _d < 3; ++_d)

        XYZGrid srcTmpXYZ(preloadgrids);
        srcTmpXYZ.open(srctempgridpath, fs::path());
        XYZGrid tgtTmpXYZ(preloadgrids);
        tgtTmpXYZ.open(tgttempgridpath, fs::path());

        Eigen::Vector3d _pair_xfm;
        Eigen::Matrix3d _pair_cov;
        bool valid = create_edge_minimize_reproject_error2(&srcTmpXYZ, &tgtTmpXYZ, poses[srcId], poses[tgtId], matches.get(), numpts, offset_n_weight.get(), invsigmasq_2, minimum_matches,
                                                           _pair_xfm.data(), _pair_cov.data(), &logger);

        double update_weight = 1. / all_cov(0, 0) + 1. / _pair_cov(0, 0);

        all_xfm = (all_xfm * 1. / all_cov(0, 0) + _pair_xfm * 1. / _pair_cov(0, 0)) / update_weight;
        all_cov = Eigen::Matrix3d::Identity() * 1. / update_weight;
      }  // for (auto& p : xyzgrid_pairs)

      memcpy(e.xfm, all_xfm.data(), sizeof(double) * 3);
      memcpy(e.cov, all_cov.data(), sizeof(double) * 9);

      double initXfm[3] = {g.nodes[tgtId].xfm[0] - g.nodes[srcId].xfm[0], g.nodes[tgtId].xfm[1] - g.nodes[srcId].xfm[1], g.nodes[tgtId].xfm[2] - g.nodes[srcId].xfm[2]};

      double diffXfm[3] = {std::abs(initXfm[0] - e.xfm[0]), std::abs(initXfm[1] - e.xfm[1]), std::abs(initXfm[2] - e.xfm[2])};

      logger.debug("Init: {} {} {}", initXfm[0], initXfm[1], initXfm[2]);
      logger.debug("Edge: {} {} {}", e.xfm[0], e.xfm[1], e.xfm[2]);
      logger.debug("Diff: {} {} {} cm {}, {}", diffXfm[0] * 100, diffXfm[1] * 100, diffXfm[2] * 100, srcName, tgtName);
//#pragma omp critical
      g.insertEdge(e);
    }
  }


  total_stat_reprojection_error.print(logger);

  ofstream ofs(outfilepath.string());
  ofs << &g;
  ofs.close();
  logger.info("Done");

  return 0;
}