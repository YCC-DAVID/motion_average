#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;
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
  fs.read((char*)datapointerout, sizeof(float) * numpt * 4);
  fs.close();

  return true;
}

namespace fs = std::filesystem;

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "exe input bin file" << endl;
    return 0;
  }
  fs::path filename(argv[1]);

  fs::path outputfile = filename;
  outputfile = outputfile.replace_extension(".txt");

  std::cout << "Input:  " << filename << endl;
  std::cout << "Output: " << outputfile << endl;

  unsigned int numpmt = 0;
  read_correspondence_pts_pix_binary_num(filename.string().c_str(), numpmt);
  std::vector<float> memdata;
  memdata.resize(numpmt * 4);

  read_correspondence_pts_pix_binary(filename.string().c_str(), memdata.data());

  std::fstream fs(outputfile, std::ios::out);
  fs.precision(15);
  fs << std::fixed;

  for (int p = 0; p < numpmt; p++) {
    fs << memdata.data()[p * 4] << " " << memdata.data()[p * 4 + 1] << " " << memdata.data()[p * 4 + 2] << " " << memdata.data()[p * 4 + 3] << std::endl;
  }
  fs.close();

  return 0;
}