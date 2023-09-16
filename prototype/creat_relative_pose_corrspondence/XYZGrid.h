#ifndef XYZGRID_H
#define XYZGRID_H

#include <gdal.h>
#include <vector>

#include <filesystem>

class XYZGrid {
 public:
  XYZGrid(bool _preload);
  ~XYZGrid();

  bool open(const std::filesystem::path gridpath[3], const std::filesystem::path raynumpath);
  bool sample(const float px, const float py, float& sx, float& sy, float& sz, uint16_t& nray) const;

  GDALDatasetH ds[4];
  GDALRasterBandH bd[4];
  bool preload = false;
  bool inited = false;
  bool hasRayNum = false;
  int width = 0, height = 0;
  std::vector<float> grids[3];
  std::vector<uint16_t> nraygrid;
};

#endif // XYZGRID_H