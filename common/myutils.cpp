#include "myutils.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "mystructs.h"

using namespace std;

bool isImg(fs::path p) {
  std::string ext = p.extension().string();
  if (ext == ".tif" || ext == ".TIF")
    return true;
  else
    return false;
}

bool isParamFile(fs::path p) {
  std::string ext = p.extension().string();
  if (ext != ".txt" && ext != ".TXT") return false;

  motionavg::Affine2D::AffineModel model;
  bool ret = model.load(p.string());
  return ret;
}

std::vector<std::string> split(const std::string& s, char seperator) {
  std::vector<std::string> output;

  std::string::size_type prev_pos = 0, pos = 0;

  while ((pos = s.find(seperator, pos)) != std::string::npos) {
    std::string substring(s.substr(prev_pos, pos - prev_pos));

    output.push_back(substring);

    prev_pos = ++pos;
  }

  output.push_back(s.substr(prev_pos, pos - prev_pos));  // Last word

  return output;
}

std::vector<std::string> split2(std::string s, std::string delimiter) {
  size_t pos_start = 0, pos_end, delim_len = delimiter.length();
  std::string token;
  std::vector<std::string> res;

  while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos) {
    token = s.substr(pos_start, pos_end - pos_start);
    pos_start = pos_end + delim_len;
    res.push_back(token);
  }

  res.push_back(s.substr(pos_start));
  return res;
}
