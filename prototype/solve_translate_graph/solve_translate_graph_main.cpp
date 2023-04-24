#include <spdlog/fmt/ostr.h>
#include <spdlog/spdlog.h>

#include <Eigen/Dense>
#include <algorithm>
#include <cxxopts.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "MoAve_0.h"
#include "MoAve_Affine.h"
#include "mystructs.h"
#include "myutils.h"

using namespace std;
namespace logg = spdlog;
using namespace motionavg::TranslateND;
using GraphT = TranslateGraph3WithGCP;

constexpr double GCPWeightFactor = 0.05;
GraphT solve_translate_MoAve_0_direct(const GraphT& g) {
  GraphT og = g;

  std::vector<std::vector<double>> NodePara;
  std::vector<std::vector<ConG_Unit>> Con_Graph;
  std::vector<std::vector<UnaryG_Unit>> Unary_Graph;

  size_t n_node = g.nodes.size();
  size_t n_pairedge = g.edges.size();
  size_t n_gcpedge = g.gcplinks.size();
  bool hasGCP = n_gcpedge > 0;

  NodePara.resize(n_node);
  for (size_t i = 0; i < n_node; ++i) {
    const GraphT::Node& n = g.nodes[i];
    NodePara[i].resize(3);
    std::copy_n(n.xfm, 3, NodePara[i].data());
  }

  Con_Graph.resize(n_node);
  Unary_Graph.resize(n_node);

  // Add pairwise constraints
  for (size_t i = 0; i < n_pairedge; ++i) {
    ConG_Unit u;
    const GraphT::Edge& e = g.edges[i];
    u.RefID = e.target;
    u.OtherID = e.source;
    u.RelParaVec.resize(3);
    std::copy_n(e.xfm, 3, u.RelParaVec.data());

    Eigen::Matrix<double, 3, 3> covMat(e.cov);
    Eigen::Matrix<double, 3, 3> PMat = covMat.inverse();

    Eigen::JacobiSVD<Eigen::Matrix<double, 3, 3>> svd(PMat, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Vector<double, 3> s_sqrt = svd.singularValues().cwiseSqrt();
    Eigen::Matrix<double, 3, 3> PMatSqrt = svd.matrixU() * s_sqrt.asDiagonal() * svd.matrixV().transpose();

    u.InvCovMatrix.resize(9);
    std::copy_n(PMatSqrt.data(), 9, u.InvCovMatrix.data());

    Con_Graph[e.target].push_back(u);
  }

  // Add direct constraints to node
  for (size_t i = 0; i < n_gcpedge; ++i) {
    UnaryG_Unit u;
    const GraphT::GCPLink& l = g.gcplinks[i];
    const GraphT::GCPPoint& gcppt = g.gcps[l.gcpid];
    u.RefID = l.viewid;
    u.AbsParaVec.resize(3);
    u.AbsParaVec[0] = gcppt.x - l.dx;
    u.AbsParaVec[1] = gcppt.y - l.dy;
    u.AbsParaVec[2] = gcppt.z - l.dz;

    u.InvCovMatrix.resize(9);

    Eigen::Matrix<double, 3, 3> stdMat;
    stdMat.setZero();
    stdMat(0, 0) = sqrt(std::pow(gcppt.ex, 2) + std::pow(l.ex, 2));
    stdMat(1, 1) = sqrt(std::pow(gcppt.ey, 2) + std::pow(l.ey, 2));
    stdMat(2, 2) = sqrt(std::pow(gcppt.ez, 2) + std::pow(l.ez, 2));
    Eigen::Matrix<double, 3, 3> PMatSqrt = stdMat.inverse() * GCPWeightFactor;
    std::copy_n(PMatSqrt.data(), 9, u.InvCovMatrix.data());

    Unary_Graph[u.RefID].push_back(u);
  }

  std::vector<int> RefNode;
  if (!hasGCP) RefNode.push_back(0); // if GCP not appear, fix first frame to prevent rank deficiency

  MoAve_0 MA;
  MA.Initialize(NodePara, 3, Con_Graph);
  MA.SetUnaryGraph(Unary_Graph);
  MA.Set_Reference_Nodes(RefNode);
  MA.Set_Max_Iter_Num(200);
  MA.Set_ConvergenceEpslon(1e-14);
  MA.DirectSolver_Run();
  // Write out

  for (size_t i = 0; i < n_node; ++i) {
    GraphT::Node& n = og.nodes[i];
    std::copy_n(MA._Node_Para_list[i].data(), 3, n.xfm);
  }
  return og;
}

GraphT solve_translate_MoAve_0_iter(const GraphT& g) {
  GraphT og = g;

  std::vector<std::vector<double>> NodePara;
  std::vector<std::vector<ConG_Unit>> Con_Graph;

  size_t n_node = g.nodes.size();
  size_t n_edge = g.edges.size();

  NodePara.resize(n_node);
  for (size_t i = 0; i < n_node; ++i) {
    const GraphT::Node& n = g.nodes[i];
    NodePara[i].resize(3);
    std::copy_n(n.xfm, 3, NodePara[i].data());
  }

  Con_Graph.resize(n_node);

  for (size_t i = 0; i < n_edge; ++i) {
    ConG_Unit u;
    const GraphT::Edge& e = g.edges[i];
    u.RefID = e.target;
    u.OtherID = e.source;
    u.RelParaVec.resize(3);
    std::copy_n(e.xfm, 3, u.RelParaVec.data());

    Eigen::Matrix<double, 3, 3> covMat(e.cov);
    Eigen::Matrix<double, 3, 3> PMat = covMat.inverse();

    Eigen::JacobiSVD<Eigen::Matrix<double, 3, 3>> svd(PMat, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Vector<double, 3> s_sqrt = svd.singularValues().cwiseSqrt();
    Eigen::Matrix<double, 3, 3> PMatSqrt = svd.matrixU() * s_sqrt.asDiagonal() * svd.matrixV().transpose();

    u.InvCovMatrix.resize(9);
    std::copy_n(PMatSqrt.data(), 9, u.InvCovMatrix.data());

    Con_Graph[e.target].push_back(u);
  }

  std::vector<int> RefNode;
  RefNode.push_back(0);

  MoAve_0 MA;
  MA.Initialize(NodePara, 3, Con_Graph);
  MA.Set_Reference_Nodes(RefNode);
  MA.Set_ConvergenceEpslon(1e-7);
  MA.Set_Max_Iter_Num(100);
  MA.IterSolver_Run();
  // Write out

  for (size_t i = 0; i < n_node; ++i) {
    GraphT::Node& n = og.nodes[i];
    std::copy_n(MA._Node_Para_list[i].data(), 3, n.xfm);
  }
  return og;
}

void evaluate(GraphT& g, bool verbose) {
  double pairwise_error[3] = {0};
  auto logger=spdlog::get("LOG");
  logger->info("==== Evaluating Pairwise ====");
  for (int ei = 0; ei < g.edges.size(); ++ei) {
    GraphT::Edge& e = g.edges[ei];
    GraphT::Node& srcN = g.nodes[e.source];
    GraphT::Node& tgtN = g.nodes[e.target];

    double initXfm[3] = {tgtN.xfm[0] - srcN.xfm[0], tgtN.xfm[1] - srcN.xfm[1], tgtN.xfm[2] - srcN.xfm[2]};

    double diffXfm[3] = {std::abs(initXfm[0] - e.xfm[0]), std::abs(initXfm[1] - e.xfm[1]), std::abs(initXfm[2] - e.xfm[2])};
    for (int i = 0; i < 3; ++i) pairwise_error[i] += diffXfm[i];
    if (verbose) {
      string sname = fs::path(srcN.name).stem().string();
      string tname = fs::path(tgtN.name).stem().string();
      logger->info("src {} tgt {}", sname, tname);
      logger->debug("Init: {} {} {}", initXfm[0], initXfm[1], initXfm[2]);
      logger->debug("Edge: {} {} {}", e.xfm[0], e.xfm[1], e.xfm[2]);
      logger->info("Diff: {} {} {} cm", diffXfm[0] * 100, diffXfm[1] * 100, diffXfm[2] * 100);
    }
  }
  for (int i = 0; i < 3; ++i) pairwise_error[i] /= g.edges.size();
  logger->info("Mean Error: {} {} {} cm", pairwise_error[0] * 100, pairwise_error[1] * 100, pairwise_error[2] * 100);

  logger->info("==== Evaluating GCP ====");
  double gcp_error[3] = {0};
  for (int ei=0;ei<g.gcplinks.size();++ei) {
    GraphT::GCPLink& l = g.gcplinks[ei];
    GraphT::Node& n = g.nodes[l.viewid];
    GraphT::GCPPoint& gcppt = g.gcps[l.gcpid];

    double initGCP[3] = {n.xfm[0] + l.dx, n.xfm[1] + l.dy, n.xfm[2] + l.dz};
    double diffGCP[3] = {std::abs(initGCP[0] - gcppt.x), std::abs(initGCP[1] - gcppt.y), std::abs(initGCP[2] - gcppt.z)};
    for (int i = 0; i < 3; ++i) gcp_error[i] += diffGCP[i];
    if (verbose) {
      string gcpname = gcppt.name;
      string viewname = fs::path(n.name).stem().string();
      logger->info("View {} and GCP {}", viewname, gcpname);
      logger->debug("Init: {} {} {}", initGCP[0], initGCP[1], initGCP[2]);
      logger->debug("Edge: {} {} {}", gcppt.x, gcppt.y, gcppt.z);
      logger->info("Diff: {} {} {} cm", diffGCP[0] * 100, diffGCP[1] * 100, diffGCP[2] * 100);
    }
  }

  for (int i = 0; i < 3; ++i) gcp_error[i] /= g.gcplinks.size();
  logger->info("Mean Error: {} {} {} cm", gcp_error[0] * 100, gcp_error[1] * 100, gcp_error[2] * 100);
}

void graph_normalize_pairwise_cov(GraphT& g) {
  auto logger = spdlog::get("LOG");
  std::vector<double> covs;
  covs.reserve(g.edges.size());
  for (size_t i = 0; i < g.edges.size(); ++i) covs.push_back(g.edges[i].cov[0]);
  auto beg = covs.begin();
  auto end = covs.end();
  auto mid = beg + g.edges.size() / 2;
  std::nth_element(beg, mid, end);

  auto median = *mid;
  logger->critical("Median Cov {}", median);
  for (size_t i = 0; i < covs.size(); ++i)
    for (int _d = 0; _d < 3; ++_d) g.edges[i].cov[_d * 4] /= median;
}

cxxopts::Options parseOptions(std::string exepath = "") {
  std::string exename = fs::path(exepath).filename().string();
  cxxopts::Options options(exename, "create_graph from folder");
  // clang-format off
  options.add_options()
    ("i,input_graph", "input graph", cxxopts::value<std::string>())
	("o,output_name", "output filename", cxxopts::value<std::string>()->default_value("graph.txt"))
	("h,help", "Print help")
	("v,verbose", "verbose level (trace - 0, debug - 1, info - 2, warn - 3, error - 4, critical - 5, off - 6)", cxxopts::value<int>()->default_value("2"))
	;
  // clang-format on

  options.parse_positional({"input_graph", "output_name"});
  options.positional_help("input_graph output_name [options]");
  return options;
}

#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>

int main(int argc, char** argv) {
  cxxopts::Options options = parseOptions(argv[0]);
  cxxopts::ParseResult args = options.parse(argc, argv);

  if (args.count("help") != 0) {
    cout << options.help() << endl;
    return 0;
  }

  if (args.count("input_graph") == 0) {
    cout << options.help() << endl;
    return 0;
  }

  fs::path input_graph(args["input_graph"].as<string>());
  fs::path output_name(args["output_name"].as<string>());
  if (output_name.is_relative()) output_name = input_graph.parent_path() / output_name;
  
  auto console_sink = std::make_shared<spdlog::sinks::stderr_color_sink_mt>();
  console_sink->set_level(spdlog::level::info);

  fs::path log_path = output_name.parent_path() / fmt::format("log_solve_{}.log", output_name.stem().string());
  auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(log_path.string());
  file_sink->set_level(spdlog::level::trace);

  std::vector<spdlog::sink_ptr> sinks;
  sinks.push_back(console_sink);
  sinks.push_back(file_sink);

  auto logger = std::make_shared<spdlog::logger>("LOG", begin(sinks), end(sinks));
  spdlog::register_logger(logger);
  logger->set_level(static_cast<logg::level::level_enum>(args["verbose"].as<int>()));

  logger->info("Input:  {}", input_graph.string());
  logger->info("Output: {}", output_name.string());

  if (!fs::exists(input_graph)) {
    logger->error("File not found: {}", input_graph.string());
    cout << options.help() << endl;
    return 0;
  }

  ifstream ifs(input_graph.string());

  GraphT g;
  ifs >> &g;
  graph_normalize_pairwise_cov(g);

  evaluate(g, false);
  g.rebase(0);
  GraphT direct = solve_translate_MoAve_0_direct(g);
  direct.rebase(-1);
  ofstream ofs(output_name.string());
  ofs << &direct;
  ofs.close();
  evaluate(direct, false);

  // GraphT iter = solve_translate_MoAve_0_iter(g);
  // iter.rebase(-1);
  // ofs.open(itergraphpath);
  // ofs << iter;
  // ofs.close();

  logger->info("Done");

  return 0;
}