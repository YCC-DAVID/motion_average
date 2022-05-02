#include <omp.h>
#include <spdlog/spdlog.h>

#include <Eigen/Dense>
#include <cxxopts.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <pdal/PointLayout.hpp>
#include <pdal/PointTable.hpp>
#include <pdal/PointView.hpp>
#include <pdal/io/BufferReader.hpp>
#include <pdal/pdal.hpp>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "../create_translate_graph_MSP/XYZGrid.h"
#include "gdal.h"
#include "mystructs.h"
#include "myutils.h"
#include "qin_io.h"

namespace fs = std::filesystem;
namespace logg = spdlog;
using namespace std;
using namespace motionavg::TranslateND;
using GraphT = TranslateGraph<3>;

cxxopts::Options parseOptions(std::string exepath = "") {
  std::string exename = fs::path(exepath).filename().string();
  cxxopts::Options options(exename, "apply translate graph from folder");
  // clang-format off
  options.add_options()
    ("i,input_folder", "input folder", cxxopts::value<std::string>())
	("graphfile", "translate_graph_file", cxxopts::value<std::string>())
	("o,output_folder", "output folder", cxxopts::value<std::string>()->default_value(""))
	("inqinv2", "input qinv2file", cxxopts::value<std::string>())
	("outqinv2", "output qinv2file", cxxopts::value<std::string>()->default_value(""))
	("las", "output las", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
	("interval", "interval of sampling", cxxopts::value<int>()->default_value("16"))
	("minrays", "minimum num of rays", cxxopts::value<int>()->default_value("3"))
	("h,help", "Print help")
	("v,verbose", "verbose level (trace - 0, debug - 1, info - 2, warn - 3, error - 4, critical - 5, off - 6)", cxxopts::value<int>()->default_value("2"))
	;
  // clang-format on

  options.parse_positional({"input_folder", "graphfile", "output_folder"});
  options.positional_help("input_folder graphfile output_folder [options]");
  return options;
}

vector<QinPose> readQinv2File(std::string input_qin) {
  ifstream ifs(input_qin);
  size_t num_pose;
  ifs >> num_pose;
  vector<QinPose> poses(num_pose);
  for (size_t i = 0; i < num_pose; ++i) ifs >> poses[i];
  return poses;
}

void writeQinv2File(std::string output_qin, const vector<QinPose>& poses) {
  ofstream ofs(output_qin);
  ofs << std::fixed << std::setprecision(15);
  size_t num_pose = poses.size();
  ofs << num_pose << endl;
  for (size_t i = 0; i < num_pose; ++i) ofs << poses[i];
}

int main(int argc, char** argv) {
  using namespace pdal;

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
  fs::path graphfilepath(args["graphfile"].as<string>());
  fs::path outdirpath(args["output_folder"].as<string>());
  int interval = args["interval"].as<int>();
  int minimum_rays = args["minrays"].as<int>();
  bool write_las = args["las"].as<bool>();
  if (graphfilepath.is_relative()) graphfilepath = rootdirpath / graphfilepath;

  if (outdirpath.empty()) outdirpath = graphfilepath.parent_path() / graphfilepath.filename().stem();
  if (!fs::exists(outdirpath)) fs::create_directories(outdirpath);

  GraphT g;
  ifstream ifs(graphfilepath.string());
  ifs >> g;
  ifs.close();

  if (args.count("inqinv2") > 0) {
    fs::path inqin(args["inqinv2"].as<string>());
    fs::path outqin(args["outqinv2"].as<string>());
    if (outqin.empty()) outqin = outdirpath / inqin.filename();
    if (outqin.is_relative()) outqin = outdirpath / outqin;

    if (fs::exists(inqin)) {
      vector<QinPose> poses = readQinv2File(inqin.string());
      if (poses.size() != g.nodes.size()) {
        spdlog::error("graph node and qinfile mismatch!");
      } else {
        if (!fs::exists(outqin.parent_path())) fs::create_directories(outqin.parent_path());
        for (size_t i = 0; i < poses.size(); ++i) {
          poses[i].x = g.nodes[i].xfm[0];
          poses[i].y = g.nodes[i].xfm[1];
          poses[i].z = g.nodes[i].xfm[2];
        }
        writeQinv2File(outqin.string(), poses);
        spdlog::info("Write new QinFile: {}", outqin.string());
      }
    } else {
      spdlog::error("Qinv2 not found: {}", inqin.string());
    }
  } else {
    spdlog::info("Skip process qin file");
  }

  if (!write_las) return 0;
  //// Genreate Las file
  fs::path tempfolderdirpath = rootdirpath / "Point_clouds" / "temp_folder_cluster";

#pragma omp parellel for
  for (int ni = 0; ni < g.nodes.size(); ++ni) {
    string imgname = fs::path(g.nodes[ni].name).stem().string();  // if (srcName != "024_059_id2571c81607_121250_Right") continue;
    fs::path gridpath[3] = {tempfolderdirpath / imgname / fmt::format("{}_Xgrid.tif", imgname), tempfolderdirpath / imgname / fmt::format("{}_Ygrid.tif", imgname),
                            tempfolderdirpath / imgname / fmt::format("{}_Zgrid.tif", imgname)};
    fs::path raynumpath = tempfolderdirpath / imgname / fmt::format("{}_ray_num.tif", imgname);
    for (int _d = 0; _d < 3; ++_d)
      if (!fs::exists(gridpath[_d])) {
        spdlog::error("File Not Found: {}", gridpath[_d].string());
        return 1;
      }
    if (!fs::exists(raynumpath)) {
      spdlog::error("File Not Found: {}", raynumpath.string());
      return 1;
    }
    string outlaspath = (outdirpath / (imgname + "_point_cloud.las")).string();

    PointTable table;
    table.layout()->registerDim(Dimension::Id::X);
    table.layout()->registerDim(Dimension::Id::Y);
    table.layout()->registerDim(Dimension::Id::Z);

    PointViewPtr view(new PointView(table));

    XYZGrid grid(false);
    grid.open(gridpath, raynumpath);

    for (int _x = 0, cnt = 0; _x < grid.width; _x += interval)
      for (int _y = 0; _y < grid.height; _y += interval) {
        float px = _x;
        float py = _y;
        float lX, lY, lZ;
        uint16_t nray;
        if (!grid.sample(px, py, lX, lY, lZ, nray)) continue;
        if (nray < minimum_rays) continue;
        double gX, gY, gZ;
        gX = g.nodes[ni].xfm[0] + lX;
        gY = g.nodes[ni].xfm[1] + lY;
        gZ = g.nodes[ni].xfm[2] + lZ;

        view->setField(pdal::Dimension::Id::X, cnt, gX);
        view->setField(pdal::Dimension::Id::Y, cnt, gY);
        view->setField(pdal::Dimension::Id::Z, cnt, gZ);
        ++cnt;
      }

    pdal::Options options;
    options.add("filename", outlaspath);
    pdal::BufferReader reader;
    reader.addView(view);

    StageFactory factory;

    // Set second argument to 'true' to let factory take ownership of
    // stage and facilitate clean up.
    Stage* writer = factory.createStage("writers.las");

    writer->setInput(reader);
    writer->setOptions(options);
    writer->prepare(table);
    writer->execute(table);

    spdlog::info("Output {}", outlaspath);
  }

  return 0;
}