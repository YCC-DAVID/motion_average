#include "myutils.h"
#include "mystructs.h"
#include <iostream>
#include <iomanip>
#include <cassert>
#include <fstream>
#include <spdlog/spdlog.h>

namespace fs = std::filesystem;
namespace logg = spdlog;
using namespace std;
using namespace motionavg::Affine2D;



void evaluate_posegraph(const BundleGraph& gt, PoseGraph& pg)
{
	double sse = 0;
	double count = 0;
	for (const auto& e : gt.edges)
	{
		assert(gt.nodes[e.source].name == pg.nodes[e.source].name);
		assert(gt.nodes[e.target].name == pg.nodes[e.target].name);
		
		auto& ns = pg.nodes[e.source];
		auto& nt = pg.nodes[e.target];
		
		double edge_sse = 0;
		double tmpS[2], tmpT[2];
		for (const auto& tpair : e.tiepoints)
		{
			XfmApply(ns.poseXfm, tpair.source[0], tpair.source[1], tmpS[0], tmpS[1]);
			XfmApply(nt.poseXfm, tpair.target[0], tpair.target[1], tmpT[0], tmpT[1]);
			edge_sse += std::pow(tmpT[0] - tmpS[0], 2) + std::pow(tmpT[1] - tmpS[1], 2);
		}
		sse += edge_sse;
		count += e.tiepoints.size();
	}
	double rmse = std::sqrt(sse / count);
	cout << "Point RMSE " << rmse <<" m" << endl;
}

void evaluate_bundlegraph(const BundleGraph& gt, BundleGraph& pg)
{
	double sse = 0;
	double count = 0;
	for (const auto& e : gt.edges)
	{
		assert(gt.nodes[e.source].name == pg.nodes[e.source].name);
		assert(gt.nodes[e.target].name == pg.nodes[e.target].name);

		auto& ns = pg.nodes[e.source];
		auto& nt = pg.nodes[e.target];

		double edge_sse = 0;
		double tmpS[2], tmpT[2];
		for (const auto& tpair : e.tiepoints)
		{
			XfmApply(ns.poseXfm, tpair.source[0], tpair.source[1], tmpS[0], tmpS[1]);
			XfmApply(nt.poseXfm, tpair.target[0], tpair.target[1], tmpT[0], tmpT[1]);
			edge_sse += std::pow(tmpT[0] - tmpS[0], 2) + std::pow(tmpT[1] - tmpS[1], 2);
		}
		sse += edge_sse;
		count += e.tiepoints.size();
	}
	double rmse = std::sqrt(sse / count);
	cout << "Point RMSE " << rmse << " m" << endl;
}


void evaluate(const BundleGraph& gt, string graphfile)
{
	ifstream ifs(graphfile);
	string identifier;
	std::getline(ifs, identifier);
	ifs.seekg(0);
	if (identifier == "PoseGraph")
	{
		PoseGraph pg;
		ifs >> pg;
		cout << "Evaluating Pose Graph " << graphfile << endl;
		evaluate_posegraph(gt, pg);
	}
	else if (identifier == "BundleGraph")
	{
		BundleGraph bg;
		ifs >> bg;
		cout << "Evaluating Bundle Graph " << graphfile << endl;
		evaluate_bundlegraph(gt, bg);
	}
	else {
		cerr << "Unknown format: " << identifier << endl;
	}
}

int main(int argc, char** argv)
{
	string gtfile = argv[1];
	vector<string> solfiles;
	for (int i = 2; i < argc; ++i)
		if(fs::exists(argv[i])) solfiles.push_back(argv[i]);
	if (!fs::exists(argv[1])) return 2;
	if (solfiles.empty()) return 2;

	BundleGraph gtgraph;
	ifstream ifs(gtfile);
	ifs >> gtgraph;
	ifs.close();
	
	cout << "Ground Truth: " << gtfile << endl;
	cout << "--------------\n";
	for (size_t i = 0; i < solfiles.size(); ++i)
	{
		evaluate(gtgraph, solfiles[i]);
	}
	return 0;
}