#include <gdal.h>
#include <cpl_string.h>
#include <filesystem>
#include <iostream>
#include <fstream>

#include "myutils.h"
#include "mystructs.h"


using namespace std;
namespace fs = std::filesystem;
using namespace motionavg::Affine2D;

bool ensureDir(fs::path p)
{
	if (!fs::exists(p))
		return fs::create_directories(p);
	else
		return true;
}

int main(int argc, char** argv)
{
	GDALAllRegister();
	string graphfilepath = argv[1];
	string outdirpath = argv[2];


	XfmGraph* graph;

	ifstream ifs(graphfilepath);
	string identifier;
	std::getline(ifs, identifier);
	ifs.seekg(0);
	if (identifier == "PoseGraph")
	{
		auto* g = new PoseGraph();
		ifs >> *g;
		graph = g;
	}
	else if (identifier == "BundleGraph") {
		auto* g = new BundleGraph();
		ifs >> *g;
		graph = g;
	}
	else {
		cerr << "Unknown format: " << identifier << endl;
		return 1;
	}

	fs::path rootdir(graph->basepath);
	fs::path outdir(outdirpath);
	if (!fs::exists(outdir))
		fs::create_directories(outdir);
	size_t num_nodes = graph->nodes.size();

	GDALDriverH driver = GDALGetDriverByName("GTiff");
	char** papszCreateOptions = NULL;
	//papszCreateOptions = CSLSetNameValue(papszCreateOptions, "TFW", "YES");
	for (size_t i = 0; i < num_nodes; ++i)
	{
		XfmGraph::Node& n = graph->nodes[i];
		fs::path inputfile = (rootdir / fs::path(n.path));
		fs::path outputfile = (outdir / fs::path(n.path));

		double globalXfm[6];
		XfmComposite(graph->baseGeo, n.poseXfm, globalXfm);

		if (!ensureDir(outputfile.parent_path()))
		{
			cerr << "Cannot create folder for " << outputfile << endl;
		}

		if (fs::exists(outputfile))
		{
			cout << "update " << outputfile << endl;
			auto dstDS = GDALOpen(outputfile.string().c_str(), GA_Update);
			GDALSetGeoTransform(dstDS, globalXfm);
			GDALClose(dstDS);
		}
		else {
			cout << inputfile << " -> " << outputfile << endl;
			auto srcDS = GDALOpen(inputfile.string().c_str(), GA_ReadOnly);
			auto dstDS = GDALCreateCopy(driver, outputfile.string().c_str(), srcDS, FALSE, papszCreateOptions, NULL, NULL);
			GDALSetGeoTransform(dstDS, globalXfm);
			GDALClose(dstDS);
			GDALClose(srcDS);
		}
	}

	delete graph;

	return 0;
}