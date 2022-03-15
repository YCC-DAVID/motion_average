#include "myutils.h"
#include "mystructs.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

using namespace std;

bool isImg(fs::path p)
{
	std::string ext = p.extension().string();
	if (ext == ".tif" || ext == ".TIF")
		return true;
	else
		return false;
}


bool isParamFile(fs::path p)
{
	std::string ext = p.extension().string();
	if (ext != ".txt" && ext != ".TXT")
		return false;

	motionavg::AffineModel model;
	bool ret = model.load(p.string());
	return ret;
}

std::vector<std::string> split(const std::string& s, char seperator)
{
	std::vector<std::string> output;

	std::string::size_type prev_pos = 0, pos = 0;

	while ((pos = s.find(seperator, pos)) != std::string::npos)
	{
		std::string substring(s.substr(prev_pos, pos - prev_pos));

		output.push_back(substring);

		prev_pos = ++pos;
	}

	output.push_back(s.substr(prev_pos, pos - prev_pos)); // Last word

	return output;
}
