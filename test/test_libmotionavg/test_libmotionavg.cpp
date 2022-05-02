#include <vector>
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include "MoAve_0.h"
#include "MoAve_Affine.h"
#include "mystructs.h"
#include "../test_common.h"

struct MoAve_0_Syn100_BadInit_TestsFixture {

	MoAve_0_Syn100_BadInit_TestsFixture() {
		generate_0order_graph(NodePara, Con_Graph, num_nodes, num_para_per_node, num_neighbor_radius);
		DistortNodePara.assign(NodePara.begin(), NodePara.end());
		DistortCon_Graph.assign(Con_Graph.begin(), Con_Graph.end());
		RefNode.push_back(0);
		distort_0order_graph(DistortNodePara, DistortCon_Graph, RefNode, 1000, 0.1);
	}
	int num_nodes = 100;
	int num_para_per_node = 3;
	int num_neighbor_radius = 6;
	std::vector<int> RefNode;
	std::vector<std::vector<double>> NodePara, DistortNodePara;
	std::vector<std::vector<ConG_Unit>> Con_Graph, DistortCon_Graph;
};


TEST_CASE_FIXTURE(MoAve_0_Syn100_BadInit_TestsFixture, "DirectSolver wo edge noise")
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
TEST_CASE_FIXTURE(MoAve_0_Syn100_BadInit_TestsFixture, "DirectSolver w edge noise"){
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
TEST_CASE_FIXTURE(MoAve_0_Syn100_BadInit_TestsFixture, "IterSolver wo edge noise")
{
	MoAve_0 M0;
	M0.Initialize(DistortNodePara, num_para_per_node, Con_Graph);
	M0.Set_Max_Iter_Num(10000);
	M0.Set_Reference_Nodes(RefNode);
	M0.Set_ConvergenceEpslon(1e-17);
	bool ret = M0.IterSolver_Run();

	REQUIRE(ret);
	for (int ni = 0; ni < num_nodes; ++ni)
		for (int j = 0; j < num_para_per_node; ++j)
			CHECK(NodePara[ni][j] == doctest::Approx(M0._Node_Para_list[ni][j]).epsilon(1e-2));
}

TEST_CASE_FIXTURE(MoAve_0_Syn100_BadInit_TestsFixture, "IterSolver w edge noise")
{
	MoAve_0 M0;
	M0.Initialize(DistortNodePara, num_para_per_node, DistortCon_Graph);
	M0.Set_Max_Iter_Num(10000);
	M0.Set_Reference_Nodes(RefNode);
	M0.Set_ConvergenceEpslon(1e-17);
	bool ret = M0.IterSolver_Run();

	REQUIRE(ret);
	for (int ni = 0; ni < num_nodes; ++ni)
		for (int j = 0; j < num_para_per_node; ++j)
			CHECK(NodePara[ni][j] == doctest::Approx(M0._Node_Para_list[ni][j]).epsilon(0.5));
}


struct MoAve_0_Syn100_GoodInit_TestsFixture {

	MoAve_0_Syn100_GoodInit_TestsFixture() {
		generate_0order_graph(NodePara, Con_Graph, num_nodes, num_para_per_node, num_neighbor_radius);
		DistortNodePara.assign(NodePara.begin(), NodePara.end());
		DistortCon_Graph.assign(Con_Graph.begin(), Con_Graph.end());
		RefNode.push_back(0);
		distort_0order_graph(DistortNodePara, DistortCon_Graph, RefNode, 1, 0.1);
	}
	int num_nodes = 100;
	int num_para_per_node = 3;
	int num_neighbor_radius = 6;
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
	M0.Set_Max_Iter_Num(10000);
	M0.Set_Reference_Nodes(RefNode);
	M0.Set_ConvergenceEpslon(1e-17);
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
	M0.Set_Max_Iter_Num(10000);
	M0.Set_Reference_Nodes(RefNode);
	M0.Set_ConvergenceEpslon(1e-17);
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
}

struct MoAve_Affine_Syn100_BadInit_TestsFixture {

	MoAve_Affine_Syn100_BadInit_TestsFixture() {
		generate_affine_graph(NodePara, Con_Graph, num_nodes, num_neighbor_radius);
		DistortNodePara.assign(NodePara.begin(), NodePara.end());
		DistortCon_Graph.assign(Con_Graph.begin(), Con_Graph.end());
		RefNode.push_back(0);
		distort_affine_graph(DistortNodePara, DistortCon_Graph, RefNode, 0.1, 1000, 1e-4, 0.1);
	}
	int num_nodes = 100;
	int num_para_per_node = 6;
	int num_neighbor_radius = 6;
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

TEST_CASE_FIXTURE(MoAve_Affine_Syn100_BadInit_TestsFixture, "DirectSolver w edge noise") {
	MoAve_Affine M0;
	M0.Initialize(DistortNodePara, num_para_per_node, DistortCon_Graph);
	M0.Set_Max_Iter_Num(100);
	M0.Set_Reference_Nodes(RefNode);
	bool ret = M0.DirectSolver_Run();

	REQUIRE(ret);
	for (int ni = 0; ni < num_nodes; ++ni)
		for (int j = 0; j < num_para_per_node; ++j)
			CHECK(NodePara[ni][j] == doctest::Approx(M0._Node_Para_list[ni][j]).epsilon(0.2));
}


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


struct MoAve_0_2DSyn100_GoodInit_TestsFixture {

	MoAve_0_2DSyn100_GoodInit_TestsFixture() {
		generate_0order_graph(NodePara, Con_Graph, num_nodes, num_para_per_node, num_neighbor_radius);
		DistortNodePara.assign(NodePara.begin(), NodePara.end());
		DistortCon_Graph.assign(Con_Graph.begin(), Con_Graph.end());
		RefNode.push_back(0);
		distort_0order_graph(DistortNodePara, DistortCon_Graph, RefNode, 5, 0.1);
	}
	int num_nodes = 20;
	int num_para_per_node = 2;
	int num_neighbor_radius = 6;
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
	M0.Initialize(DistortNodePara, num_para_per_node, Con_Graph);
	M0.Set_Max_Iter_Num(10000);
	M0.Set_Reference_Nodes(RefNode);
	M0.Set_ConvergenceEpslon(1e-17);
	bool ret = M0.IterSolver_Run();

	REQUIRE(ret);
	for (int ni = 0; ni < num_nodes; ++ni)
		for (int j = 0; j < num_para_per_node; ++j)
			CHECK(NodePara[ni][j] == doctest::Approx(M0._Node_Para_list[ni][j]).epsilon(0.05));
}
