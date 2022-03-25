#include <spdlog/spdlog.h>
#include <spdlog/fmt/ostr.h>
#include <cxxopts.hpp>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "mystructs.h"
#include "myutils.h"


#include "MoAve_0.h"
#include "MoAve_Affine.h"

#include <Eigen/Dense>

using namespace std;
using namespace motionavg;


PoseGraph solve_posegraph_MoAve_Affine_Iter(const PoseGraph& g)
{
	PoseGraph og = g;

	std::vector<std::vector<double>> NodePara;
	std::vector<std::vector<ConG_Unit>> Con_Graph;

	size_t n_node = g.nodes.size();
	size_t n_edge = g.edges.size();

	NodePara.resize(n_node);
	for (size_t i = 0; i < n_node; ++i)
	{
		const PoseGraph::Node& n = g.nodes[i];
		NodePara[i].resize(6);
		std::copy_n(n.poseXfm, 6, NodePara[i].data());
	}

	Con_Graph.resize(n_node);

	for (size_t i = 0; i < n_edge; ++i)
	{
		ConG_Unit u;
		const PoseGraph::Edge& e = g.edges[i];
		u.RefID = e.target;
		u.OtherID = e.source;
		u.RelParaVec.resize(6);
		std::copy_n(e.regXfm, 6, u.RelParaVec.data());

		Eigen::Matrix<double, 6, 6> covMat(e.covXfm);
		Eigen::Matrix<double, 6, 6> PMat = covMat.inverse();

		Eigen::JacobiSVD< Eigen::Matrix<double, 6, 6>> svd(PMat, Eigen::ComputeFullU | Eigen::ComputeFullV);
		Eigen::Vector<double, 6> s_sqrt = svd.singularValues().cwiseSqrt();
		Eigen::Matrix<double, 6, 6> PMatSqrt = svd.matrixU() * s_sqrt.asDiagonal() * svd.matrixV().transpose();

		u.InvCovMatrix.resize(36);
		std::copy_n(PMatSqrt.data(), 36, u.InvCovMatrix.data());

		Con_Graph[e.target].push_back(u);
	}

	
	std::vector<int> RefNode;
	RefNode.push_back(3);

	MoAve_Affine MA;
	MA.Initialize(NodePara, 6, Con_Graph);
	MA.Set_Reference_Nodes(RefNode);
	MA.Set_Max_Iter_Num(50);
	MA.IterSolver_Run();
	//MA.DirectSolver_Run();
	// Write out

	for (size_t i = 0; i < n_node; ++i)
	{
		PoseGraph::Node& n = og.nodes[i];
		std::copy_n(MA._Node_Para_list[i].data(), 6, n.poseXfm);
	}
	return og;
}

int main(int argc, char** argv)
{
	string ingraph = "N:\\motion_avg_whole\\graph.txt";
	string outgraph = "N:\\motion_avg_whole\\graph_analytic.txt";
	ifstream ifs(ingraph);
	string typeline;
	getline(ifs, typeline);
	ifs.seekg(0);
	if (typeline == "PoseGraph")
	{
		PoseGraph g;
		ifs >> g;
		spdlog::info("PoseGraph {} {} {}", g.basepath, g.nodes.size(), g.edges.size());
		ifs.close();

		PoseGraph og = solve_posegraph_MoAve_Affine_Iter(g);

		ofstream ofs(outgraph);
		ofs << og;
		ofs.close();
	}
	else if (typeline == "BundleGraph")
	{
		BundleGraph g;
		ifs >> g;
		spdlog::info("BundleGraph {} {} {}", g.basepath, g.nodes.size(), g.edges.size());
		ifs.close();
	}
	else {
		spdlog::error("Unknown Graph Type: {}", typeline);
		ifs.close();
		return 1; 
	}
	

	return 0;
}