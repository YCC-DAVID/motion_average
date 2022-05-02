#include <spdlog/spdlog.h>

#include <cassert>
#include <cxxopts.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "mystructs.h"
#include "myutils.h"

namespace fs = std::filesystem;
namespace logg = spdlog;
using namespace std;
using namespace motionavg::Affine2D;

void evaluate(const BundleGraph& gt, XfmGraph& g) {
  double sse = 0;
  double count = 0;
  for (const auto& e : gt.edges) {
    assert(gt.nodes[e.source].name == g.nodes[e.source].name);
    assert(gt.nodes[e.target].name == g.nodes[e.target].name);

    auto& ns = g.nodes[e.source];
    auto& nt = g.nodes[e.target];

    double edge_sse = 0;
    double tmpS[2], tmpT[2];
    for (const auto& tpair : e.tiepoints) {
      XfmApply(ns.poseXfm, tpair.source[0], tpair.source[1], tmpS[0], tmpS[1]);
      XfmApply(nt.poseXfm, tpair.target[0], tpair.target[1], tmpT[0], tmpT[1]);
      edge_sse += std::pow(tmpT[0] - tmpS[0], 2) + std::pow(tmpT[1] - tmpS[1], 2);
    }
    sse += edge_sse;
    count += e.tiepoints.size();
  }
  double rmse = std::sqrt(sse / count);
  logg::info("Point RMSE {} m", rmse);
}

// ------------------------------------------
enum GRAPH_TYPE { POSE, BUNDLE, CASTPOSE, NUM_TYPE };

cxxopts::Options parseOptions(std::string exepath = "") {
  std::string exename = fs::path(exepath).filename().string();
  cxxopts::Options options(exename, "evaluate solution quality with raw tie point pair. Compute RMSE in world reference system.");
  // clang-format off
  options.add_options()
    ("gt", "input BundleGraph file as Ground Truth", cxxopts::value<std::string>())
	("graphlist", "input graph files for evaluation.", cxxopts::value<std::vector<std::string>>())
	("h,help", "Print help")
	("v,verbose", "verbose level (trace - 0, debug - 1, info - 2, warn - 3, error - 4, critical - 5, off - 6)", cxxopts::value<int>()->default_value("2"))
	;
  // clang-format on

  options.parse_positional({"gt", "graphlist"});
  options.positional_help("gt solution1,solution2,solution3,solution4,...");
  return options;
}

int main(int argc, char** argv) {
  cxxopts::Options options = parseOptions(argv[0]);
  cxxopts::ParseResult args = options.parse(argc, argv);
  spdlog::set_level(static_cast<logg::level::level_enum>(args["verbose"].as<int>()));

  if (args.count("help") != 0) {
    cout << options.help() << endl;
    ;
    return 0;
  }

  if (args.count("gt") == 0) {
    logg::error("GroundTruth not found.");
    cout << options.help() << endl;
    return 1;
  } else if (!fs::exists(args["gt"].as<string>())) {
    logg::error("GroundTruth not exists");
    cout << options.help() << endl;
    return 1;
  }

  if (args.count("graphlist") == 0) {
    logg::error("graphlist not found.");
    cout << options.help() << endl;
    return 1;
  }

  string gtfile = args["gt"].as<string>();
  vector<string> inputgraphlist = args["graphlist"].as<vector<string>>();

  for (auto p : args.unmatched()) {
    inputgraphlist.push_back(p);
  }

  vector<string> graphlist;

  for (auto p : inputgraphlist) {
    if (fs::exists(p))
      graphlist.push_back(p);
    else
      logg::warn("Skip {}: not exists", p);
  }

  BundleGraph gtgraph;
  ifstream ifs(gtfile);
  ifs >> gtgraph;
  ifs.close();

  logg::info("Ground Truth: {}", gtfile);
  logg::info("--------------");

  for (auto p : graphlist) {
    XfmGraph g;
    ifs.open(p);
    ifs >> g;
    ifs.close();
    logg::info("{}:", p);
    evaluate(gtgraph, g);
  }
  return 0;
}