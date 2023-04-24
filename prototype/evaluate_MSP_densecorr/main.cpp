#include <spdlog/spdlog.h>
#include <cxxopts.hpp>

#include <string>
#include <vector>
#include <filesystem>
#include <fstream>
#include "qin_io.h"
#include "XYZGrid.h"
#include <Eigen/Eigen>

namespace fs = std::filesystem;
using namespace std;


template <typename T>
class StatVar {
  std::string _name;
  T _E;
  T _sqE;
  uint32_t _cnt;

 public:
  StatVar(std::string name) : _name(name) { _cnt = 0; }
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

  void print(std::ostream& out) { out << "E(" << _name << ")=" << getE().transpose() << ", D(" << _name << ")=" << getVar().transpose() << " of " << getCount() << " measurements" << std::endl; }
};


bool read_correspondence_pts_pix_binary_num(std::string filename, unsigned int& numpts) {
  std::fstream fs(filename, std::ios::in | std::ios::binary);
  if (!fs.is_open()) {
    std::cout << "file opening error" << std::endl;
    return false;
  }
  fs.read((char*)&numpts, sizeof(unsigned int));
  fs.close();
}
//
bool read_correspondence_pts_pix_binary(std::string filename, float* datapointerout) {
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

vector<QinPose> readQinv2File(std::string input_qin) {
  ifstream ifs(input_qin);
  size_t num_pose;
  ifs >> num_pose;
  vector<QinPose> poses(num_pose);
  for (size_t i = 0; i < num_pose; ++i) ifs >> poses[i];
  return poses;
}

int main(int argc, char** argv) { 
    GDALAllRegister();

    fs::path qin_path(R"(M:\MoAve\Dublin\Area1\images_undist\pose.qinv2)");
    fs::path res_path(R"(M:\MoAve\Dublin\Area1\res)");
    fs::path Neighborhood_file_path = res_path / "Point_clouds" / "point_clouds" / "Neighborhood_file.bin";
    fs::path temp_folder_path = res_path / "Point_clouds" / "temp_folder_cluster";

    auto poses = readQinv2File(qin_path.string());

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

    std::vector<std::vector<int>> decoded_neighbors(numimg);
    for (int i = 0; i < numimg; ++i)  // i-th image
    {
    const int* ids = &neighborbin[i * (maxpair + 1)];
    for (const int* _pj = ids + 1; _pj < ids + maxpair + 1; ++_pj) {
        if (*_pj == -1) break;
        decoded_neighbors[i].push_back(*_pj);
    }
    }


    for (int srcId=0;srcId<numimg;++srcId) {
      for (int tgtId : decoded_neighbors[srcId]) {
        std::string srcName = fs::path(poses[srcId].imgname).stem().string();
        std::string tgtName = fs::path(poses[tgtId].imgname).stem().string();
        fs::path corr_bin_path = temp_folder_path / srcName / fmt::format("{}_{}_correspondence_pix.bin", srcName, tgtName);
        if(!fs::exists(corr_bin_path)) {
          spdlog::error("CORR.BIN NOT FOUND");
          continue;
        }
        spdlog::info(corr_bin_path.stem().string());

        
        unsigned int numpts;
        read_correspondence_pts_pix_binary_num(corr_bin_path.string(), numpts);
        spdlog::info("numpts {}", numpts);
        std::unique_ptr<float> matches(new float[7 * numpts]);
        read_correspondence_pts_pix_binary(corr_bin_path.string(), matches.get());


        // clang-format off
        fs::path tmpgridpath[3] = {
            temp_folder_path / srcName / fmt::format("{}_{}_temgridX.tif", srcName, tgtName),
            temp_folder_path / srcName / fmt::format("{}_{}_temgridY.tif", srcName, tgtName),
            temp_folder_path / srcName / fmt::format("{}_{}_temgridZ.tif", srcName, tgtName)
        };
        // clang-format on
        XYZGrid tmpXYZ(false);
        tmpXYZ.open(tmpgridpath, fs::path());


        auto srcK = poses[srcId].GetK();
        auto srcR = poses[srcId].GetR();
        auto srcC = poses[srcId].GetC();

        auto tgtK = poses[tgtId].GetK();
        auto tgtR = poses[tgtId].GetR();
        auto tgtC = poses[tgtId].GetC();

        for (int ptId = 0; ptId<numpts;++ptId) {
        
            const float* basePtr = matches.get() + ptId * 7;
            Eigen::Vector2d srcPt{basePtr[0], basePtr[1]};
            Eigen::Vector2d tgtPt{basePtr[2], basePtr[3]};
            Eigen::Vector3d srcObj{basePtr[4], basePtr[5], basePtr[6]};

            float _gridX, _gridY, _gridZ;
            uint16_t dummy;
            tmpXYZ.sample(srcPt[0], srcPt[1], _gridX, _gridY, _gridZ, dummy);


            Eigen::Vector3d tgtObj = srcObj + srcC - tgtC;
            Eigen::Vector3d tgtCam = tgtK * tgtR * tgtObj;
            tgtCam /= tgtCam.z();
            spdlog::info("[{},{}]: {} {} {} from grid {} {} {}", srcPt[0], srcPt[1], srcObj[0], srcObj[1], srcObj[2], _gridX, _gridY, _gridZ);
            spdlog::info("[{},{}]: {},{} from 3d {} {}",srcPt[0], srcPt[1], tgtPt[0], tgtPt[1], tgtCam[0], tgtCam[1]);
        }

        /*
        // Now, make combination of src-neighbors and tgt-neighbors, except src-tgt
        std::vector<std::pair<int, int>> xyzgrid_pairs;
        for (int adj_srcId : decoded_neighbors[srcId])
          for (int adj_tgtId : decoded_neighbors[tgtId]) {
            if (adj_srcId != tgtId && adj_tgtId != srcId) xyzgrid_pairs.emplace_back(adj_srcId, adj_tgtId);
          }
          */



        break;
      }
      break;
    }








    
    return 0;
}