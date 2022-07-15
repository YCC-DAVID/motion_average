#include "XYZGrid.h"



#include <gdal.h>
#include <spdlog/spdlog.h>

#include <filesystem>

namespace fs = std::filesystem;

XYZGrid::XYZGrid(bool _preload) : preload(_preload) {}

XYZGrid::~XYZGrid() {
  for (int _d = 0; _d < 4; ++_d) GDALClose(ds[_d]);
}

bool XYZGrid::open(const fs::path gridpath[3], const fs::path raynumpath) {
  for (int _d = 0; _d < 3; ++_d) {
    ds[_d] = GDALOpen(gridpath[_d].string().c_str(), GA_ReadOnly);
    bd[_d] = GDALGetRasterBand(ds[_d], 1);
    if (GDALGetRasterDataType(bd[_d]) != GDT_Float32) {
      spdlog::error("Data Type is Not Float32: {}", gridpath[_d].string());
      return false;
    }
  }
  if (!raynumpath.empty()) {
    ds[3] = GDALOpen(raynumpath.string().c_str(), GA_ReadOnly);
    bd[3] = GDALGetRasterBand(ds[3], 1);
    if (GDALGetRasterDataType(bd[3]) != GDT_UInt16) {
      spdlog::error("Data Type is Not UInt16: {}", raynumpath.string());
      return false;
    }
    hasRayNum = true;
  }

  width = GDALGetRasterXSize(ds[0]);
  height = GDALGetRasterYSize(ds[0]);

  if (preload) {
    size_t bufsize_ = width * height;
    for (int _d = 0; _d < 3; ++_d) {
      grids[_d].resize(bufsize_);
      GDALRasterIO(bd[_d], GF_Read, 0, 0, width, height, grids[_d].data(), width, height, GDT_Float32, sizeof(float), width * sizeof(float));
    }
    if (hasRayNum) {
      nraygrid.resize(bufsize_);
      GDALRasterIO(bd[3], GF_Read, 0, 0, width, height, nraygrid.data(), width, height, GDT_UInt16, sizeof(uint16_t), width * sizeof(uint16_t));
    }
  }

  inited = true;
  return inited;
}

  bool XYZGrid::sample(const float& px, const float& py, float& sx, float& sy, float& sz, uint16_t& nray) const {
  float* res[3] = {&sx, &sy, &sz};
  bool valid = true;
  if (preload) {
    int c0 = px;
    int c1 = px + 1;
    int r0 = py;
    int r1 = py + 1;
    float cw = px - c0;
    float rw = py - r0;
    for (int _d = 0; _d < 3; ++_d) {
      // v00 ---------- v01
      //  |			   |
      //  |    .(px, py) |
      //  |			   |
      // v10 ---------- v11
      const float& v00 = grids[_d][r0 * width + c0];
      const float& v01 = grids[_d][r0 * width + c1];
      const float& v10 = grids[_d][r1 * width + c0];
      const float& v11 = grids[_d][r1 * width + c1];
      float vr0 = v00 * (1 - cw) + v01 * cw;
      float vr1 = v10 * (1 - cw) + v01 * cw;
      *res[_d] = vr0 * (1 - rw) + vr1 * rw;
      if (std::isnan(*res[_d]) || std::isinf(*res[_d])) valid = false;
    }
    {
      const uint16_t& v00 = nraygrid[r0 * width + c0];
      const uint16_t& v01 = nraygrid[r0 * width + c1];
      const uint16_t& v10 = nraygrid[r1 * width + c0];
      const uint16_t& v11 = nraygrid[r1 * width + c1];
      if (hasRayNum)
        nray = std::min<uint16_t>(v00, std::min<uint16_t>(v01, std::min<uint16_t>(v10, v11)));
      else
        nray = 0;
    }
  } else {
    GDALRasterIOExtraArg args;
    args.bFloatingPointWindowValidity = true;
    args.dfXOff = px;
    args.dfYOff = py;
    args.dfXSize = 1;
    args.dfYSize = 1;
    args.eResampleAlg = GRIORA_Bilinear;
    for (int _d = 0; _d < 3; ++_d) {
      GDALRasterIO(bd[_d], GF_Read, px, py, 1, 1, res[_d], 1, 1, GDT_Float32, sizeof(float), sizeof(float));
      if (std::isnan(*res[_d]) || std::isinf(*res[_d])) valid = false;
    }
    if (hasRayNum)
      GDALRasterIO(bd[3], GF_Read, px, py, 1, 1, &nray, 1, 1, GDT_UInt16, sizeof(uint16_t), sizeof(uint16_t));
    else
      nray = 0;
  }
  return valid;
}