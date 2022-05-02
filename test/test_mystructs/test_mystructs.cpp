#include "mystructs.h"
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include <spdlog/spdlog.h>

using namespace motionavg::TranslateND;
using namespace motionavg::Affine2D;

void genreate_random_pose_graph(PoseGraph& g)
{
	g.basepath = std::to_string(rand());
	int num_nodes = 100;
	int num_neighbor_radius = 6;
	g.nodes.resize(num_nodes);
	for (int ni = 0; ni < num_nodes; ++ni)
	{
		g.nodes[ni].name = std::to_string(rand());
		g.nodes[ni].path = std::to_string(rand());
		g.nodes[ni].poseXfm[0] = 0.01 * (rand() % 10000);
		g.nodes[ni].poseXfm[3] = 0.01 * (rand() % 10000);
		g.nodes[ni].poseXfm[1] = 1.0 + 0.2 * float(rand()) / float(RAND_MAX);
		g.nodes[ni].poseXfm[2] = 0.0 + 0.1 * float(rand()) / float(RAND_MAX);
		g.nodes[ni].poseXfm[4] = 0.0 + 0.1 * float(rand()) / float(RAND_MAX);
		g.nodes[ni].poseXfm[5] = 1.0 + 0.2 * float(rand()) / float(RAND_MAX);
		g.nodes[ni].left = rand();
		g.nodes[ni].right = rand();
		g.nodes[ni].top = rand();
		g.nodes[ni].bottom = rand();
	}
	g.edges.clear();
	for (int ni = 0; ni < num_nodes; ++ni)
	{
		for (int j = -num_neighbor_radius; j <= num_neighbor_radius; ++j)
		{
			if (j == 0) continue;
			PoseGraph::Edge edge;
			edge.name = std::to_string(rand());
			edge.target = ni;
			edge.source = (ni + num_nodes - j) % num_nodes;

			for (int l = 0; l < 6; ++l)
			{
				edge.covXfm[6 * l + l] = 0.1*(rand()%100);
			}
			double tem[6];
			XfmInv(g.nodes[edge.target].poseXfm, tem);
			XfmComposite(tem, g.nodes[edge.source].poseXfm, edge.regXfm);
			g.edges.push_back(edge);
		}
	}
}

template<int N>
void generate_random_translate_graph(TranslateGraph<N>& g)
{
	using GraphT = TranslateGraph<N>;
	g.baseID = -1;
	g.basepath = std::to_string(rand());
	for (int i = 0; i < N; ++i)
		g.baseGeo[i] = 0.01 * (rand() % 1000);
	
	for (int i = 0; i < 100; ++i)
	{
		typename GraphT::Node n;
		n.name = std::to_string(rand());
		n.path = std::to_string(rand());
		for (int j = 0; j < GraphT::DIM; ++j)
		{
			n.xfm[j] = 0.1 * double(rand() % 1000);
		}
		g.insertNode(n);
	}

	for (int i = 0; i < 200; ++i)
	{
		typename GraphT::Edge e;
		e.name = std::to_string(rand());
		e.source = rand() % 10;
		e.target = rand() % 10;
		for (int j = 0; j < GraphT::DIM; ++j)
		{
			e.xfm[j] = 0.1 * double(rand() % 1000);
		}for (int j = 0; j < GraphT::DIM * GraphT::DIM; ++j)
		{
			e.cov[j] = 0.1 * double(rand() % 1000);
		}
		g.insertEdge(e);
	}
}

TEST_CASE("Translate Graph 2D")
{
	char tmpfilename[256];
	tmpnam_s(tmpfilename);
	spdlog::info("Temp file: {}", tmpfilename);

	std::ifstream ifs;
	std::ofstream ofs;

	using GraphT = TranslateGraph<2>;
	GraphT g, g2;
	generate_random_translate_graph(g);
	
	ofs.open(tmpfilename);
	ofs << g;
	ofs.close();

	g.rebase(2);
	g.rebase(-1);
	g.rebase(1);


	ifs.open(tmpfilename);
	ifs >> g2;
	ifs.close();
	g2.rebase(-1);
	g2.rebase(5);
	g2.rebase(-1);
	g2.rebase(1);

	CHECK(g.baseID == g2.baseID);
	for (int i = 0; i < GraphT::DIM; ++i)
		CHECK(g.baseGeo[i] == doctest::Approx(g2.baseGeo[i]));
	CHECK(g.basepath == g2.basepath);
	REQUIRE(g.nodes.size()==g2.nodes.size());
	REQUIRE(g.edges.size() == g2.edges.size());
	for (int i = 0; i < g.nodes.size(); ++i)
	{
		CHECK(g.nodes[i].name == g2.nodes[i].name);
		CHECK(g.nodes[i].path == g2.nodes[i].path);
		for (int j = 0; j < GraphT::DIM; ++j)
			CHECK(g.nodes[i].xfm[j] == doctest::Approx(g2.nodes[i].xfm[j]).epsilon(1e-7));
	}
	for (int i = 0; i < g.edges.size(); ++i)
	{
		CHECK(g.edges[i].name == g2.edges[i].name);
		for (int j = 0; j < GraphT::DIM; ++j)
			CHECK(g.edges[i].xfm[j] == doctest::Approx(g2.edges[i].xfm[j]).epsilon(1e-7));
		for (int j = 0; j < GraphT::DIM* GraphT::DIM; ++j)
			CHECK(g.edges[i].cov[j] == doctest::Approx(g2.edges[i].cov[j]).epsilon(1e-7));
	}
	std::remove(tmpfilename);
}


TEST_CASE("Translate Graph 3D")
{
	char tmpfilename[256];
	tmpnam_s(tmpfilename);
	spdlog::info("Temp file: {}", tmpfilename);

	std::ifstream ifs;
	std::ofstream ofs;

	using GraphT = TranslateGraph<3>;
	GraphT g, g2;
	generate_random_translate_graph(g);

	ofs.open(tmpfilename);
	ofs << g;
	ofs.close();

	g.rebase(2);
	g.rebase(1);


	ifs.open(tmpfilename);
	ifs >> g2;
	ifs.close();
	g2.rebase(5);
	g2.rebase(1);

	CHECK(g.baseID == g2.baseID);
	for (int i = 0; i < GraphT::DIM; ++i)
		CHECK(g.baseGeo[i] == doctest::Approx(g2.baseGeo[i]));
	CHECK(g.basepath == g2.basepath);
	REQUIRE(g.nodes.size() == g2.nodes.size());
	REQUIRE(g.edges.size() == g2.edges.size());
	for (int i = 0; i < g.nodes.size(); ++i)
	{
		CHECK(g.nodes[i].name == g2.nodes[i].name);
		CHECK(g.nodes[i].path == g2.nodes[i].path);
		for (int j = 0; j < GraphT::DIM; ++j)
			CHECK(g.nodes[i].xfm[j] == doctest::Approx(g2.nodes[i].xfm[j]).epsilon(1e-7));
	}
	for (int i = 0; i < g.edges.size(); ++i)
	{
		CHECK(g.edges[i].name == g2.edges[i].name);
		for (int j = 0; j < GraphT::DIM; ++j)
			CHECK(g.edges[i].xfm[j] == doctest::Approx(g2.edges[i].xfm[j]).epsilon(1e-7));
		for (int j = 0; j < GraphT::DIM * GraphT::DIM; ++j)
			CHECK(g.edges[i].cov[j] == doctest::Approx(g2.edges[i].cov[j]).epsilon(1e-7));
	}
	std::remove(tmpfilename);
}


TEST_CASE("Translate Graph")
{
	char tmpfilename[256];
	tmpnam_s(tmpfilename);
	spdlog::info("Temp file: {}", tmpfilename);

	std::ifstream ifs;
	std::ofstream ofs;

	using GraphT = PoseGraph;
	GraphT g, g2;
	genreate_random_pose_graph(g);
	

	ofs.open(tmpfilename);
	ofs << g;
	ofs.close();

	g.rebase(3);
	g.rebase(-1);
	g.rebase(6);

	ifs.open(tmpfilename);
	ifs >> g2;
	ifs.close();
	g2.rebase(-1);
	g2.rebase(9);
	g2.rebase(-1);
	g2.rebase(1);
	g2.rebase(6);

	CHECK(g.baseID == g2.baseID);
	CHECK(g.basepath == g2.basepath);
	for (int i = 0; i < GraphT::DIM; ++i)
		CHECK(g.baseGeo[i] == doctest::Approx(g2.baseGeo[i]));
	REQUIRE(g.nodes.size() == g2.nodes.size());
	REQUIRE(g.edges.size() == g2.edges.size());
	for (int i = 0; i < g.nodes.size(); ++i)
	{
		CHECK(g.nodes[i].name == g2.nodes[i].name);
		CHECK(g.nodes[i].path == g2.nodes[i].path);
		for (int j = 0; j < GraphT::DIM; ++j)
			CHECK(g.nodes[i].poseXfm[j] == doctest::Approx(g2.nodes[i].poseXfm[j]).epsilon(1e-7));
	}
	for (int i = 0; i < g.edges.size(); ++i)
	{
		CHECK(g.edges[i].name == g2.edges[i].name);
		for (int j = 0; j < GraphT::DIM; ++j)
			CHECK(g.edges[i].regXfm[j] == doctest::Approx(g2.edges[i].regXfm[j]).epsilon(1e-7));
		for (int j = 0; j < GraphT::DIM * GraphT::DIM; ++j)
			CHECK(g.edges[i].covXfm[j] == doctest::Approx(g2.edges[i].covXfm[j]).epsilon(1e-7));
	}
	std::remove(tmpfilename);
}