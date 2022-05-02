#include <spdlog/spdlog.h>

#include <Eigen/Dense>
#include <cxxopts.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "gdal.h"
#include "mystructs.h"
#include "myutils.h"

namespace fs = std::filesystem;
namespace logg = spdlog;
using namespace std;
using namespace motionavg::Affine2D;
// ------------------------------------------

bool search_qins_files(fs::path rootdir, std::unordered_map<string, fs::path>& imgindex, std::unordered_map<string, PairwiseData>& paramindex) {
  // --------------------------------
  // Collect necessary data from Dr. Qin's output
  logg::info("Searching in directory ... [{0}]", rootdir.string());

  std::vector<fs::path> imglist;
  std::vector<fs::path> paramlist;

  for (const auto& f : fs::recursive_directory_iterator(rootdir)) {
    if (f.is_directory()) continue;
    if (isImg(f))
      imglist.push_back(f);
    else if (isParamFile(f))
      paramlist.push_back(f);
  }

  logg::info("Found {0} Images and {1} parameters", imglist.size(), paramlist.size());

  for (const auto& p : imglist) imgindex[p.filename().stem().string()] = p;
  logg::info("Img name index: {0}", imgindex.size());
  if (imgindex.size() != imglist.size()) logg::warn("Warning: images with same name detected!");

  for (const auto& p : paramlist) {
    PairwiseData model;
    model.load(p.string());
    paramindex[model.affine.filename] = model;
  }
  logg::info("Param name index: {0}", paramindex.size());
  if (paramindex.size() != paramindex.size()) logg::warn("Warning: transformations with same name detected!");
  return true;
}
PoseGraph build_pose_graph(fs::path rootdir) {
  std::unordered_map<string, fs::path> imgindex;
  std::unordered_map<string, PairwiseData> paramindex;
  search_qins_files(rootdir, imgindex, paramindex);
  // --------------------------------------------
  // Build Pose Graph
  logg::info("Building PoseGraph ...");

  PoseGraph pgraph;
  pgraph.basepath = rootdir.string();

  std::unordered_map<string, size_t> imgname_to_nid;

  for (const auto& kv : paramindex) {
    string name = kv.second.affine.filename;
    string target_img = kv.second.affine.target_file;
    string source_img = kv.second.affine.source_file;

    int target_nid = -1, source_nid = -1;
    if (imgname_to_nid.find(target_img) != imgname_to_nid.end()) {
      target_nid = imgname_to_nid[target_img];
    } else {
      PoseGraph::Node n;
      n.name = target_img;
      n.path = fs::relative(imgindex[target_img], rootdir).string();
      GDALDatasetH pds = GDALOpen(imgindex[target_img].string().c_str(), GA_ReadOnly);
      n.left = 0;
      n.right = GDALGetRasterXSize(pds);
      n.top = 0;
      n.bottom = GDALGetRasterYSize(pds);
      double pGeoTrans[6];
      GDALGetGeoTransform(pds, pGeoTrans);
      std::copy_n(pGeoTrans, 6, n.poseXfm);
      GDALClose(pds);
      target_nid = pgraph.insertNode(n);
      imgname_to_nid[target_img] = target_nid;
    }

    if (imgname_to_nid.find(source_img) != imgname_to_nid.end()) {
      source_nid = imgname_to_nid[source_img];
    } else {
      PoseGraph::Node n;
      n.name = source_img;
      n.path = fs::relative(imgindex[source_img], rootdir).string();
      GDALDatasetH pds = GDALOpen(imgindex[source_img].string().c_str(), GA_ReadOnly);
      n.left = 0;
      n.right = GDALGetRasterXSize(pds);
      n.top = 0;
      n.bottom = GDALGetRasterYSize(pds);
      double pGeoTrans[6];
      GDALGetGeoTransform(pds, pGeoTrans);
      std::copy_n(pGeoTrans, 6, n.poseXfm);
      GDALClose(pds);

      source_nid = pgraph.insertNode(n);
      imgname_to_nid[source_img] = source_nid;
    }

    PoseGraph::Edge e;
    e.name = name;
    e.target = target_nid;
    e.source = source_nid;
    // GeoParam is a adjustment for Right Pose. The relative position should be:
    const PairwiseData& pairdata = kv.second;
    const AffineModel& affinemodel = pairdata.affine;
    double originSourceGeo[6], resampleSourceGeo[6], composeSourceGeo[6];
    std::copy_n(pgraph.nodes[source_nid].poseXfm, 6, originSourceGeo);
    std::copy_n(originSourceGeo, 6, resampleSourceGeo);
    XfmApply(affinemodel.geoParam, originSourceGeo[0] - affinemodel.pixelParam[0], originSourceGeo[3] + affinemodel.pixelParam[3], resampleSourceGeo[0], resampleSourceGeo[3]);
    XfmComposite(resampleSourceGeo, affinemodel.pixelParam, composeSourceGeo);  // This is the final affine transformation from pixel location to registered UTM

    double invTargetGeo[6];
    double invTargetSourceGeo[6];
    XfmInv(pgraph.nodes[target_nid].poseXfm, invTargetGeo);
    XfmComposite(invTargetGeo, composeSourceGeo, invTargetSourceGeo);
    std::copy_n(invTargetSourceGeo, 6, e.regXfm);

    pgraph.insertEdge(e);

    logg::debug("Pose Edge: {0}", e.name);
    logg::debug("From {0}[{1}] to {2}[{3}]", e.target, pgraph.nodes[target_nid].name, e.source, pgraph.nodes[source_nid].name);
  }
  return pgraph;
}

BundleGraph build_bundle_graph(fs::path rootdir) {
  std::unordered_map<string, fs::path> imgindex;
  std::unordered_map<string, PairwiseData> paramindex;
  search_qins_files(rootdir, imgindex, paramindex);

  logg::info("Building BundleGraph ...");

  BundleGraph bgraph;
  bgraph.basepath = rootdir.string();

  std::unordered_map<string, size_t> imgname_to_nid;

  for (const auto& kv : paramindex) {
    string name = kv.second.affine.filename;
    string target_img = kv.second.affine.target_file;
    string source_img = kv.second.affine.source_file;

    int target_nid = -1, source_nid = -1;
    if (imgname_to_nid.find(target_img) != imgname_to_nid.end()) {
      target_nid = imgname_to_nid[target_img];
    } else {
      BundleGraph::Node n;
      n.name = target_img;
      n.path = fs::relative(imgindex[target_img], rootdir).string();
      GDALDatasetH pds = GDALOpen(imgindex[target_img].string().c_str(), GA_ReadOnly);
      n.left = 0;
      n.right = GDALGetRasterXSize(pds);
      n.top = 0;
      n.bottom = GDALGetRasterYSize(pds);
      double pGeoTrans[6];
      GDALGetGeoTransform(pds, pGeoTrans);
      std::copy_n(pGeoTrans, 6, n.poseXfm);
      GDALClose(pds);
      target_nid = bgraph.insertNode(n);
      imgname_to_nid[target_img] = target_nid;
    }

    if (imgname_to_nid.find(source_img) != imgname_to_nid.end()) {
      source_nid = imgname_to_nid[source_img];
    } else {
      BundleGraph::Node n;
      n.name = source_img;
      n.path = fs::relative(imgindex[source_img], rootdir).string();
      GDALDatasetH pds = GDALOpen(imgindex[source_img].string().c_str(), GA_ReadOnly);
      n.left = 0;
      n.right = GDALGetRasterXSize(pds);
      n.top = 0;
      n.bottom = GDALGetRasterYSize(pds);
      double pGeoTrans[6];
      GDALGetGeoTransform(pds, pGeoTrans);
      std::copy_n(pGeoTrans, 6, n.poseXfm);
      GDALClose(pds);

      source_nid = bgraph.insertNode(n);
      imgname_to_nid[source_img] = source_nid;
    }

    BundleGraph::Edge e;
    e.name = name;
    e.target = target_nid;
    e.source = source_nid;
    const BundleGraph::Node& fn = bgraph.nodes[target_nid];
    const BundleGraph::Node& tn = bgraph.nodes[source_nid];
    const PairwiseData& pairdata = kv.second;
    e.tiepoints.resize(pairdata.tiepoints.size());

    std::copy(pairdata.tiepoints.begin(), pairdata.tiepoints.end(), e.tiepoints.begin());

    bgraph.insertEdge(e);

    logg::debug("Pose Edge: {0}", e.name);
    logg::debug("From {0}[{1}] to {2}[{3}]", e.target, bgraph.nodes[target_nid].name, e.source, bgraph.nodes[source_nid].name);
  }
  return bgraph;
}

PoseGraph build_pose_from_bundle_graph(const BundleGraph bgraph) {
  logg::info("Converting BundleGraph to PoseGraph ...");
  PoseGraph pgraph;
  pgraph.basepath = bgraph.basepath;
  pgraph.nodes.resize(bgraph.nodes.size());
  pgraph.edges.resize(bgraph.edges.size());
  // Cast node
  for (size_t i = 0; i < pgraph.nodes.size(); ++i) {
    const auto& bnode = bgraph.nodes[i];
    auto& pnode = pgraph.nodes[i];
    pnode.name = bnode.name;
    pnode.path = bnode.path;
    pnode.left = bnode.left;
    pnode.top = bnode.top;
    pnode.right = bnode.right;
    pnode.bottom = bnode.bottom;
    std::copy_n(bnode.poseXfm, 6, pnode.poseXfm);
  }

  vector<size_t> num_tiepoints;
  num_tiepoints.reserve(bgraph.edges.size());
  for (const auto& e : bgraph.edges) num_tiepoints.push_back(e.tiepoints.size());
  std::sort(num_tiepoints.begin(), num_tiepoints.end());
  double median_tiepoints = num_tiepoints[num_tiepoints.size() / 2];
  // Cast Edge: Estimate from tie points
  double tmpGeo[6], tmpGeo2[6];
  for (size_t i = 0; i < pgraph.edges.size(); ++i) {
    const auto& bedge = bgraph.edges[i];
    auto& pedge = pgraph.edges[i];
    pedge.name = bedge.name;
    pedge.source = bedge.source;
    pedge.target = bedge.target;

    const auto& sourcenode = pgraph.nodes[pedge.source];
    const auto& targetnode = pgraph.nodes[pedge.target];

    size_t num_pt = bedge.tiepoints.size();
    Eigen::MatrixXd A = Eigen::MatrixXd(2 * num_pt, 6), b = Eigen::MatrixXd(2 * num_pt, 1);

    A.setZero();

    XfmInv(targetnode.poseXfm, tmpGeo);
    XfmComposite(tmpGeo, sourcenode.poseXfm, tmpGeo2);

    for (size_t j = 0; j < num_pt; ++j) {
      const auto& tp = bedge.tiepoints[j];
      b(2 * j) = tp.target[0];
      b(2 * j + 1) = tp.target[1];
      A(2 * j, 0) = 1;
      A(2 * j, 1) = tp.source[0];
      A(2 * j, 2) = tp.source[1];
      A(2 * j + 1, 3) = 1;
      A(2 * j + 1, 4) = tp.source[0];
      A(2 * j + 1, 5) = tp.source[1];
    }
    Eigen::MatrixXd x = A.colPivHouseholderQr().solve(b);
    Eigen::MatrixXd cov = (A.transpose() * A / median_tiepoints).inverse();
    std::copy_n(x.data(), 6, pedge.regXfm);
    std::copy_n(cov.data(), 36, pedge.covXfm);
  }
  return pgraph;
}
// ------------------------------------------
enum GRAPH_TYPE { POSE, BUNDLE, CASTPOSE, NUM_TYPE };

cxxopts::Options parseOptions(std::string exepath = "") {
  std::string exename = fs::path(exepath).filename().string();
  cxxopts::Options options(exename, "create_graph from folder");
  // clang-format off
  options.add_options()
    ("i,input_folder", "input folder", cxxopts::value<std::string>())
    ("o,output_name", "output filename", cxxopts::value<std::string>()->default_value("graph.txt"))
    ("t,type", "type of graph, options: pose, bundle, castpose", cxxopts::value<std::string>()->default_value("castpose"))
    ("h,help", "Print help")
    ("v,verbose", "verbose level (trace - 0, debug - 1, info - 2, warn - 3, error - 4, critical - 5, off - 6)", cxxopts::value<int>()->default_value("2"))
    ;
  // clang-format on

  options.parse_positional({"input_folder", "output_name"});
  options.positional_help("input_folder output_name [options]");
  return options;
}

int main(int argc, char** argv) {
  GDALAllRegister();

  cxxopts::Options options = parseOptions(argv[0]);
  cxxopts::ParseResult args = options.parse(argc, argv);
  spdlog::set_level(static_cast<logg::level::level_enum>(args["verbose"].as<int>()));

  if (args.count("help") != 0) {
    cout << options.help() << endl;
    ;
    return 0;
  }

  fs::path rootdirpath(args["input_folder"].as<string>());
  fs::path outfilename(args["output_name"].as<string>()), outfilepath;

  if (outfilename.is_absolute())
    outfilepath = outfilename;
  else
    outfilepath = rootdirpath / outfilename;
  //

  GRAPH_TYPE type;
  string typestr = args["type"].as<string>();
  if (typestr == "pose")
    type = POSE;
  else if (typestr == "bundle")
    type = BUNDLE;
  else if (typestr == "castpose")
    type = CASTPOSE;
  else {
    type = NUM_TYPE;
    cerr << "Unrecognized graph type: " << typestr << endl;
    cerr << options.help() << endl;
    return 1;
  }

  ofstream ofs;
  switch (type) {
    case POSE: {
      PoseGraph pg = build_pose_graph(rootdirpath);
      ofs.open(outfilepath.string());
      ofs << pg;
      logg::info("Write PoseGraph to {0}.", outfilepath.string());
      break;
    }
    case BUNDLE: {
      BundleGraph bg = build_bundle_graph(rootdirpath);
      ofs.open(outfilepath.string());
      ofs << bg;
      logg::info("Write BundleGraph to {0}.", outfilepath.string());
      break;
    }
    case CASTPOSE: {
      BundleGraph bg = build_bundle_graph(rootdirpath);
      PoseGraph pg = build_pose_from_bundle_graph(bg);
      ofs.open(outfilepath.string());
      ofs << pg;
      logg::info("Write PoseGraph to {0}.", outfilepath.string());
      break;
    }
  }
  ofs.close();
  return 0;
}