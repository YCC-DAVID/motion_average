#pragma once

#include <string>
#include <vector>
#include <filesystem>

namespace fs = std::filesystem;

bool isImg(fs::path p);
bool isParamFile(fs::path p);
std::vector<std::string> split(const std::string& s, char seperator);
