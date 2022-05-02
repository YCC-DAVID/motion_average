#include <spdlog/fmt/ostr.h>
#include <spdlog/spdlog.h>

#include <Eigen/Dense>
#include <cxxopts.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include "MoAve_0.h"
#include "MoAve_Affine.h"
#include "mystructs.h"
#include "myutils.h"

using namespace std;
using namespace motionavg::TranslateND;
using GraphT = TranslateGraph<3>;

GraphT solve_translate_MoAve_0_direct(const GraphT& g) {
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
  double error[3] = {0};
  for (int ei = 0; ei < g.edges.size(); ++ei) {
    GraphT::Edge& e = g.edges[ei];
    GraphT::Node& srcN = g.nodes[e.source];
    GraphT::Node& tgtN = g.nodes[e.target];

    double initXfm[3] = {tgtN.xfm[0] - srcN.xfm[0], tgtN.xfm[1] - srcN.xfm[1], tgtN.xfm[2] - srcN.xfm[2]};

    double diffXfm[3] = {std::abs(initXfm[0] - e.xfm[0]), std::abs(initXfm[1] - e.xfm[1]), std::abs(initXfm[2] - e.xfm[2])};
    for (int i = 0; i < 3; ++i) error[i] += diffXfm[i];
    if (verbose) {
      string sname = fs::path(srcN.name).stem().string();
      string tname = fs::path(tgtN.name).stem().string();
      spdlog::info("src {} tgt {}", sname, tname);
      spdlog::debug("Init: {} {} {}", initXfm[0], initXfm[1], initXfm[2]);
      spdlog::debug("Edge: {} {} {}", e.xfm[0], e.xfm[1], e.xfm[2]);
      spdlog::info("Diff: {} {} {} cm", diffXfm[0] * 100, diffXfm[1] * 100, diffXfm[2] * 100);
    }
  }
  for (int i = 0; i < 3; ++i) error[i] /= g.edges.size();
  spdlog::info("Error: {} {} {} cm", error[0] * 100, error[1] * 100, error[2] * 100);
}

void normalize_graph(GraphT& g) {
  std::vector<double> covs;
  covs.reserve(g.edges.size());
  for (size_t i = 0; i < g.edges.size(); ++i) covs.push_back(g.edges[i].cov[0]);
  auto beg = covs.begin();
  auto end = covs.end();
  auto mid = beg + g.edges.size() / 2;
  std::nth_element(beg, mid, end);

  auto median = *mid;
  spdlog::critical("Median Cov {}", median);
  for (size_t i = 0; i < covs.size(); ++i)
    for (int _d = 0; _d < 3; ++_d) g.edges[i].cov[_d * 4] /= median;
}


int main(int argc, char** argv) {
  spdlog::set_level(spdlog::level::info);
  string ingraph = "N:\\citymapper_70\\res\\graph.txt";
  string directgraphpath = "N:\\citymapper_70\\res\\graph_direct.txt";
  string itergraphpath = "N:\\citymapper_70\\res\\graph_iter.txt";
  ifstream ifs(ingraph);

  GraphT g;
  ifs >> g;
  normalize_graph(g);

  evaluate(g, false);
  g.rebase(0);
  GraphT direct = solve_translate_MoAve_0_direct(g);
  direct.rebase(-1);
  ofstream ofs(directgraphpath);
  ofs << direct;
  ofs.close();
  evaluate(direct, false);

  // GraphT iter = solve_translate_MoAve_0_iter(g);
  // iter.rebase(-1);
  // ofs.open(itergraphpath);
  // ofs << iter;
  // ofs.close();

  spdlog::info("Done");

  return 0;
}