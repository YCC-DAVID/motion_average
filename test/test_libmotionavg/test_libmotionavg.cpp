#include <vector>
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include "MoAve_0.h"
#include "MoAve_Affine.h"
#include "mystructs.h"

using namespace motionavg;

void generate_0order_graph(std::vector<std::vector<double>>& NodePara, std::vector<std::vector<ConG_Unit>>& Con_Graph, int num_nodes = 100, int num_para_per_node=3, int num_neighbor_radius = 3)
{
	// Generate Nodes
	NodePara.resize(num_nodes);
	for (int ni = 0; ni < num_nodes; ++ni)
	{
		NodePara[ni].resize(num_para_per_node);
		for (int pi = 0; pi < num_para_per_node; ++pi)
			NodePara[ni][pi] = 0.01 * (rand() % 10000);
	}

	// Generate Edges
	Con_Graph.resize(num_nodes);
	for (int ni = 0; ni < num_nodes; ++ni)
	{
		Con_Graph[ni].clear();
		for (int j = -num_neighbor_radius; j <= num_neighbor_radius; ++j)
		{
			if (j == 0) continue;
			ConG_Unit edge;
			edge.RefID = ni;
			edge.OtherID = (ni + num_nodes - j) % num_nodes;
			edge.RelParaVec.resize(num_para_per_node);
			edge.InvCovMatrix.assign(num_para_per_node * num_para_per_node, 0);
			for (int l = 0; l < num_para_per_node; ++l)
			{
				edge.RelParaVec[l] = NodePara[edge.RefID][l] - NodePara[edge.OtherID][l];
				edge.InvCovMatrix[num_para_per_node * l + l] = 1;
			}
			
			Con_Graph[ni].push_back(edge);
		}
	}
}

void distort_0order_graph(std::vector<std::vector<double>>& NodePara, std::vector<std::vector<ConG_Unit>>& Con_Graph, const std::vector<int>& RefNode, double sigma_node, double sigma_edge)
{
	int num_nodes = NodePara.size();
	if (num_nodes == 0) return;
	int num_para_per_node = NodePara.front().size();

	if (sigma_node > 0)
	{
		for (int ni = 0; ni < num_nodes; ++ni)
		{
			if (std::find(RefNode.begin(), RefNode.end(), ni) != RefNode.end()) continue;
			for (int j = 0; j < num_para_per_node; ++j)
			{
				NodePara[ni][j] += (float)rand() / (float)RAND_MAX * sigma_node;
			}
		}
	}

	if (sigma_edge > 0)
	{
		for (int ni = 0; ni < num_nodes; ++ni)
		{
			for (int m = 0; m < Con_Graph[ni].size(); ++m)
			{
				for (int j = 0; j < num_para_per_node; ++j)
				{
					Con_Graph[ni][m].RelParaVec[j] += (float)rand() / (float)RAND_MAX * sigma_edge;
				}
			}
		}
	}
}
// --------------------------------------------------------------------------------------------
void generate_affine_graph(std::vector<std::vector<double>>& NodePara, std::vector<std::vector<ConG_Unit>>& Con_Graph, int num_nodes = 100, int num_neighbor_radius = 6)
{
	constexpr int num_para_per_node = 6;
	// Generate Nodes
	NodePara.resize(num_nodes);
	for (int ni = 0; ni < num_nodes; ++ni)
	{
		NodePara[ni].resize(num_para_per_node);
		NodePara[ni][0] = 0.01 * (rand() % 10000);
		NodePara[ni][3] = 0.01 * (rand() % 10000);
		NodePara[ni][1] = 1.0 + 0.2 * float(rand()) / float(RAND_MAX);
		NodePara[ni][2] = 0.0 + 0.1 * float(rand()) / float(RAND_MAX);
		NodePara[ni][4] = 0.0 + 0.1 * float(rand()) / float(RAND_MAX);
		NodePara[ni][5] = 1.0 + 0.2 * float(rand()) / float(RAND_MAX);
	}

	// Generate Edges
	Con_Graph.resize(num_nodes);
	for (int ni = 0; ni < num_nodes; ++ni)
	{
		Con_Graph[ni].clear();
		for (int j = -num_neighbor_radius; j <= num_neighbor_radius; ++j)
		{
			if (j == 0) continue;
			ConG_Unit edge;
			edge.RefID = ni;
			edge.OtherID = (ni + num_nodes - j) % num_nodes;
			edge.RelParaVec.resize(num_para_per_node);
			edge.InvCovMatrix.assign(num_para_per_node * num_para_per_node, 0);
			for (int l = 0; l < num_para_per_node; ++l)
			{
				edge.InvCovMatrix[num_para_per_node * l + l] = 1;
			}
			double tem[6];
			XfmInv(NodePara[edge.RefID].data(), tem);
			XfmComposite(tem, NodePara[edge.OtherID].data(), edge.RelParaVec.data());
			Con_Graph[ni].push_back(edge);
		}
	}
}

void distort_affine_graph(std::vector<std::vector<double>>& NodePara, std::vector<std::vector<ConG_Unit>>& Con_Graph, const std::vector<int>& RefNode, double sigma_node_linear, double sigma_node_translation, double sigma_edge_linear, double sigma_edge_translation)
{
	int num_nodes = NodePara.size();
	if (num_nodes == 0) return;
	int num_para_per_node = NodePara.front().size();

	if (sigma_node_linear > 0 || sigma_node_translation > 0)
	{
		for (int ni = 0; ni < num_nodes; ++ni)
		{
			if (std::find(RefNode.begin(), RefNode.end(), ni) != RefNode.end()) continue;
			
			NodePara[ni][0] += (float)rand() / (float)RAND_MAX * sigma_node_translation;
			NodePara[ni][3] += (float)rand() / (float)RAND_MAX * sigma_node_translation;

			NodePara[ni][1] += (float)rand() / (float)RAND_MAX * sigma_node_linear;
			NodePara[ni][2] += (float)rand() / (float)RAND_MAX * sigma_node_linear;
			NodePara[ni][4] += (float)rand() / (float)RAND_MAX * sigma_node_linear;
			NodePara[ni][5] += (float)rand() / (float)RAND_MAX * sigma_node_linear;
		}
	}

	if (sigma_edge_linear > 0 || sigma_edge_translation > 0)
	{
		for (int ni = 0; ni < num_nodes; ++ni)
		{
			for (int m = 0; m < Con_Graph[ni].size(); ++m)
			{
				Con_Graph[ni][m].RelParaVec[0] += (float)rand() / (float)RAND_MAX * sigma_edge_translation;
				Con_Graph[ni][m].RelParaVec[3] += (float)rand() / (float)RAND_MAX * sigma_edge_translation;

				Con_Graph[ni][m].RelParaVec[1] += (float)rand() / (float)RAND_MAX * sigma_edge_linear;
				Con_Graph[ni][m].RelParaVec[2] += (float)rand() / (float)RAND_MAX * sigma_edge_linear;
				Con_Graph[ni][m].RelParaVec[4] += (float)rand() / (float)RAND_MAX * sigma_edge_linear;
				Con_Graph[ni][m].RelParaVec[5] += (float)rand() / (float)RAND_MAX * sigma_edge_linear;
			}
		}
	}
}

// -------------------------------------------------------------------------------------------------------------------

//struct MoAve_0_Syn100_BadInit_TestsFixture {
//
//	MoAve_0_Syn100_BadInit_TestsFixture() {
//		generate_0order_graph(NodePara, Con_Graph, num_nodes, num_para_per_node);
//		DistortNodePara.assign(NodePara.begin(), NodePara.end());
//		DistortCon_Graph.assign(Con_Graph.begin(), Con_Graph.end());
//		RefNode.push_back(0);
//		distort_0order_graph(DistortNodePara, DistortCon_Graph, RefNode, 1000, 0.1);
//	}
//	int num_nodes = 100;
//	int num_para_per_node = 3;
//	std::vector<int> RefNode;
//	std::vector<std::vector<double>> NodePara, DistortNodePara;
//	std::vector<std::vector<ConG_Unit>> Con_Graph, DistortCon_Graph;
//};
//
//
//TEST_CASE_FIXTURE(MoAve_0_Syn100_BadInit_TestsFixture, "DirectSolver wo edge noise")
//{
//	MoAve_0 M0;
//	M0.Initialize(DistortNodePara, num_para_per_node, Con_Graph);
//	M0.Set_Max_Iter_Num(100);
//	M0.Set_Reference_Nodes(RefNode);
//	bool ret = M0.DirectSolver_Run();
//
//	REQUIRE(ret);
//	for (int ni = 0; ni < num_nodes; ++ni)
//		for (int j = 0; j < num_para_per_node; ++j)
//			CHECK(NodePara[ni][j] == doctest::Approx(M0._Node_Para_list[ni][j]).epsilon(1e-4));
//}
//TEST_CASE_FIXTURE(MoAve_0_Syn100_BadInit_TestsFixture, "DirectSolver w edge noise"){
//	MoAve_0 M0;
//	M0.Initialize(DistortNodePara, num_para_per_node, DistortCon_Graph);
//	M0.Set_Max_Iter_Num(100);
//	M0.Set_Reference_Nodes(RefNode);
//	bool ret = M0.DirectSolver_Run();
//
//	REQUIRE(ret);
//	for (int ni = 0; ni < num_nodes; ++ni)
//		for (int j = 0; j < num_para_per_node; ++j)
//			CHECK(NodePara[ni][j] == doctest::Approx(M0._Node_Para_list[ni][j]).epsilon(0.5));
//}
//TEST_CASE_FIXTURE(MoAve_0_Syn100_BadInit_TestsFixture, "IterSolver wo edge noise")
//{
//	MoAve_0 M0;
//	M0.Initialize(DistortNodePara, num_para_per_node, Con_Graph);
//	M0.Set_Max_Iter_Num(1000000);
//	M0.Set_Reference_Nodes(RefNode);
//	M0.Set_ConvergenceEpslon(1e-17);
//	bool ret = M0.IterSolver_Run();
//
//	REQUIRE(ret);
//	for (int ni = 0; ni < num_nodes; ++ni)
//		for (int j = 0; j < num_para_per_node; ++j)
//			CHECK(NodePara[ni][j] == doctest::Approx(M0._Node_Para_list[ni][j]).epsilon(1e-2));
//}

//TEST_CASE_FIXTURE(MoAve_0_Syn100_BadInit_TestsFixture, "IterSolver w edge noise")
//{
//	MoAve_0 M0;
//	M0.Initialize(DistortNodePara, num_para_per_node, DistortCon_Graph);
//	M0.Set_Max_Iter_Num(1000);
//	M0.Set_Reference_Nodes(RefNode);
//	M0.Set_ConvergenceEpslon(1e-17);
//	bool ret = M0.IterSolver_Run();
//
//	REQUIRE(ret);
//	for (int ni = 0; ni < num_nodes; ++ni)
//		for (int j = 0; j < num_para_per_node; ++j)
//			CHECK(NodePara[ni][j] == doctest::Approx(M0._Node_Para_list[ni][j]).epsilon(0.5));
//}

/*
struct MoAve_0_Syn100_GoodInit_TestsFixture {

	MoAve_0_Syn100_GoodInit_TestsFixture() {
		generate_0order_graph(NodePara, Con_Graph, num_nodes, num_para_per_node);
		DistortNodePara.assign(NodePara.begin(), NodePara.end());
		DistortCon_Graph.assign(Con_Graph.begin(), Con_Graph.end());
		RefNode.push_back(0);
		distort_0order_graph(DistortNodePara, DistortCon_Graph, RefNode, 1, 0.1);
	}
	int num_nodes = 100;
	int num_para_per_node = 3;
	std::vector<int> RefNode;
	std::vector<std::vector<double>> NodePara, DistortNodePara;
	std::vector<std::vector<ConG_Unit>> Con_Graph, DistortCon_Graph;

	static void write_nodes(std::string filepath, const std::vector<std::vector<double>>& nodes)
	{
		std::ofstream ofs(filepath, std::ios::trunc | std::ios::binary);
		size_t num = nodes.size();
		ofs.write(reinterpret_cast<const char*>(&num), sizeof(size_t));
		for (size_t i = 0; i < num; ++i)
			ofs.write(reinterpret_cast<const char*>(nodes[i].data()), sizeof(double) * nodes[i].size());
		ofs.close();
	}

	static void write_edges(std::string filepath, const std::vector<std::vector<ConG_Unit>>& edges)
	{
		std::ofstream ofs(filepath, std::ios::app | std::ios::binary);
		size_t num = edges.size();
		ofs.write(reinterpret_cast<const char*>(&num), sizeof(size_t));
		for (size_t i = 0; i < num; ++i)
		{
			size_t num_eg = edges[i].size();
			ofs.write(reinterpret_cast<const char*>(&num_eg), sizeof(size_t));
			for(size_t j=0; j< num_eg; ++j)
				ofs.write(reinterpret_cast<const char*>(edges[i][j].RelParaVec.data()), sizeof(double) * edges[i][j].RelParaVec.size());
		}
		ofs.close();
	}
};

TEST_CASE_FIXTURE(MoAve_0_Syn100_GoodInit_TestsFixture, "DirectSolver wo edge noise")
{
	MoAve_0 M0;
	M0.Initialize(DistortNodePara, num_para_per_node, Con_Graph);
	M0.Set_Max_Iter_Num(100);
	M0.Set_Reference_Nodes(RefNode);
	bool ret = M0.DirectSolver_Run();

	REQUIRE(ret);
	for (int ni = 0; ni < num_nodes; ++ni)
		for (int j = 0; j < num_para_per_node; ++j)
			CHECK(NodePara[ni][j] == doctest::Approx(M0._Node_Para_list[ni][j]).epsilon(1e-4));
}

TEST_CASE_FIXTURE(MoAve_0_Syn100_GoodInit_TestsFixture, "DirectSolver w edge noise")
{
	MoAve_0 M0;
	M0.Initialize(DistortNodePara, num_para_per_node, DistortCon_Graph);
	M0.Set_Max_Iter_Num(100);
	M0.Set_Reference_Nodes(RefNode);
	bool ret = M0.DirectSolver_Run();

	REQUIRE(ret);
	for (int ni = 0; ni < num_nodes; ++ni)
		for (int j = 0; j < num_para_per_node; ++j)
			CHECK(NodePara[ni][j] == doctest::Approx(M0._Node_Para_list[ni][j]).epsilon(0.5));
}

TEST_CASE_FIXTURE(MoAve_0_Syn100_GoodInit_TestsFixture, "IterSolver wo edge noise")
{
	write_nodes("N:\\gt.bin", NodePara);
	write_edges("N:\\gt.bin", Con_Graph);
	MoAve_0 M0;
	M0.Initialize(DistortNodePara, num_para_per_node, Con_Graph);
	M0.Set_Max_Iter_Num(100);
	M0.Set_Reference_Nodes(RefNode);
	bool ret = M0.IterSolver_Run();

	REQUIRE(ret);
	for (int ni = 0; ni < num_nodes; ++ni)
		for (int j = 0; j < num_para_per_node; ++j)
			CHECK(NodePara[ni][j] == doctest::Approx(M0._Node_Para_list[ni][j]).epsilon(1e-4));
}

TEST_CASE_FIXTURE(MoAve_0_Syn100_GoodInit_TestsFixture, "IterSolver w edge noise")
{
	MoAve_0 M0;
	M0.Initialize(DistortNodePara, num_para_per_node, DistortCon_Graph);
	M0.Set_Max_Iter_Num(100);
	M0.Set_Reference_Nodes(RefNode);
	bool ret = M0.IterSolver_Run();

	REQUIRE(ret);
	for (int ni = 0; ni < num_nodes; ++ni)
		for (int j = 0; j < num_para_per_node; ++j)
			CHECK(NodePara[ni][j] == doctest::Approx(M0._Node_Para_list[ni][j]).epsilon(0.5));
}

TEST_CASE("QIN_MoAve_0") {
	std::vector<std::vector<double>> NodePara;
	std::vector<std::vector<ConG_Unit>> Con_Graph;
	int num_para_per_node = 1;

	std::vector<double> curPara; curPara.push_back(0);

	NodePara.push_back(curPara);  // 1
	NodePara.push_back(curPara);  // 2
	NodePara.push_back(curPara);  // 3

	std::vector<ConG_Unit> curCGvec; curCGvec.resize(2);
	Con_Graph.push_back(curCGvec); // 1
	Con_Graph.push_back(curCGvec); // 2
	Con_Graph.push_back(curCGvec); // 3

	Con_Graph[0][0].RefID = 0;
	Con_Graph[0][0].OtherID = 1;
	Con_Graph[0][0].RelParaVec.push_back(-2);

	Con_Graph[0][1].RefID = 0;
	Con_Graph[0][1].OtherID = 2;
	Con_Graph[0][1].RelParaVec.push_back(-1);

	Con_Graph[1][0].RefID = 1;
	Con_Graph[1][0].OtherID = 0;
	Con_Graph[1][0].RelParaVec.push_back(2);

	Con_Graph[1][1].RefID = 1;
	Con_Graph[1][1].OtherID = 2;
	Con_Graph[1][1].RelParaVec.push_back(1.2);

	Con_Graph[2][0].RefID = 2;
	Con_Graph[2][0].OtherID = 0;
	Con_Graph[2][0].RelParaVec.push_back(1);

	Con_Graph[2][1].RefID = 2;
	Con_Graph[2][1].OtherID = 1;
	Con_Graph[2][1].RelParaVec.push_back(-1.2);

	std::vector<int> RefNode;
	RefNode.push_back(0);

	MoAve_0 M0;
	M0.Initialize(NodePara, num_para_per_node, Con_Graph);
	M0.Set_Reference_Nodes(RefNode);
	bool ret = M0.IterSolver_Run();

	CHECK(ret);
}*/

struct MoAve_Affine_Syn100_BadInit_TestsFixture {

	MoAve_Affine_Syn100_BadInit_TestsFixture() {
		generate_affine_graph(NodePara, Con_Graph, num_nodes);
		DistortNodePara.assign(NodePara.begin(), NodePara.end());
		DistortCon_Graph.assign(Con_Graph.begin(), Con_Graph.end());
		RefNode.push_back(0);
		distort_affine_graph(DistortNodePara, DistortCon_Graph, RefNode, 0.1, 1000, 1e-4, 0.1);
	}
	int num_nodes = 100;
	int num_para_per_node = 6;
	std::vector<int> RefNode;
	std::vector<std::vector<double>> NodePara, DistortNodePara;
	std::vector<std::vector<ConG_Unit>> Con_Graph, DistortCon_Graph;
};

TEST_CASE_FIXTURE(MoAve_Affine_Syn100_BadInit_TestsFixture, "DirectSolver wo edge noise")
{
	MoAve_Affine M0;
	M0.Initialize(DistortNodePara, num_para_per_node, Con_Graph);
	M0.Set_Max_Iter_Num(100);
	M0.Set_Reference_Nodes(RefNode);
	bool ret = M0.DirectSolver_Run();

	REQUIRE(ret);
	for (int ni = 0; ni < num_nodes; ++ni)
		for (int j = 0; j < num_para_per_node; ++j)
			CHECK(NodePara[ni][j] == doctest::Approx(M0._Node_Para_list[ni][j]).epsilon(1e-4));
}

//TEST_CASE_FIXTURE(MoAve_Affine_Syn100_BadInit_TestsFixture, "DirectSolver w edge noise") {
//	MoAve_Affine M0;
//	M0.Initialize(DistortNodePara, num_para_per_node, DistortCon_Graph);
//	M0.Set_Max_Iter_Num(100);
//	M0.Set_Reference_Nodes(RefNode);
//	bool ret = M0.DirectSolver_Run();
//
//	REQUIRE(ret);
//	for (int ni = 0; ni < num_nodes; ++ni)
//		for (int j = 0; j < num_para_per_node; ++j)
//			CHECK(NodePara[ni][j] == doctest::Approx(M0._Node_Para_list[ni][j]).epsilon(0.2));
//}


TEST_CASE_FIXTURE(MoAve_Affine_Syn100_BadInit_TestsFixture, "IterSolver wo edge noise")
{
	MoAve_Affine M0;
	M0.Initialize(DistortNodePara, num_para_per_node, Con_Graph);
	M0.Set_Max_Iter_Num(10000);
	M0.Set_Reference_Nodes(RefNode);
	M0.Set_ConvergenceEpslon(1e-25);
	bool ret = M0.IterSolver_Run();

	REQUIRE(ret);
	for (int ni = 0; ni < num_nodes; ++ni)
		for (int j = 0; j < num_para_per_node; ++j)
			CHECK(NodePara[ni][j] == doctest::Approx(M0._Node_Para_list[ni][j]).epsilon(1e-4));
}

/*
struct MoAve_0_2DSyn100_GoodInit_TestsFixture {

	MoAve_0_2DSyn100_GoodInit_TestsFixture() {
		generate_0order_graph(NodePara, Con_Graph, num_nodes, num_para_per_node, num_neighbor_radius);
		DistortNodePara.assign(NodePara.begin(), NodePara.end());
		DistortCon_Graph.assign(Con_Graph.begin(), Con_Graph.end());
		RefNode.push_back(0);
		distort_0order_graph(DistortNodePara, DistortCon_Graph, RefNode, 5, 0.1);
	}
	int num_nodes = 10;
	int num_para_per_node = 2;
	int num_neighbor_radius = 3;
	std::vector<int> RefNode;
	std::vector<std::vector<double>> NodePara, DistortNodePara;
	std::vector<std::vector<ConG_Unit>> Con_Graph, DistortCon_Graph;

	static void write_nodes(std::string filepath, const std::vector<std::vector<double>>& nodes)
	{
		std::ofstream ofs(filepath, std::ios::trunc | std::ios::binary);
		size_t num = nodes.size();
		ofs.write(reinterpret_cast<const char*>(&num), sizeof(size_t));
		for (size_t i = 0; i < num; ++i)
			ofs.write(reinterpret_cast<const char*>(nodes[i].data()), sizeof(double) * nodes[i].size());
		ofs.close();
	}

	static void write_edges(std::string filepath, const std::vector<std::vector<ConG_Unit>>& edges)
	{
		std::ofstream ofs(filepath, std::ios::app | std::ios::binary);
		size_t num = edges.size();
		ofs.write(reinterpret_cast<const char*>(&num), sizeof(size_t));
		for (size_t i = 0; i < num; ++i)
		{
			size_t num_eg = edges[i].size();
			ofs.write(reinterpret_cast<const char*>(&num_eg), sizeof(size_t));
			for (size_t j = 0; j < num_eg; ++j)
			{
				ofs.write(reinterpret_cast<const char*>(&edges[i][j].OtherID), sizeof(int));
				ofs.write(reinterpret_cast<const char*>(edges[i][j].RelParaVec.data()), sizeof(double) * edges[i][j].RelParaVec.size());
			}
		}
		ofs.close();
	}
};

TEST_CASE_FIXTURE(MoAve_0_2DSyn100_GoodInit_TestsFixture, "IterSolver wo edge noise")
{
	write_nodes("N:\\dump_graph\\gt.bin", NodePara);
	write_edges("N:\\dump_graph\\gt.bin", Con_Graph);
	MoAve_0 M0;
	M0.Initialize(DistortNodePara, num_para_per_node, DistortCon_Graph);
	M0.Set_Max_Iter_Num(200);
	M0.Set_Reference_Nodes(RefNode);
	M0.Set_ConvergenceEpslon(1e-17);
	bool ret = M0.IterSolver_Run();

	REQUIRE(ret);
	for (int ni = 0; ni < num_nodes; ++ni)
		for (int j = 0; j < num_para_per_node; ++j)
			CHECK(NodePara[ni][j] == doctest::Approx(M0._Node_Para_list[ni][j]).epsilon(0.05));
}*/
