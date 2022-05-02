#pragma once

#include "mystructs.h"

using namespace motionavg::Affine2D;

void generate_0order_graph(std::vector<std::vector<double>>& NodePara, std::vector<std::vector<ConG_Unit>>& Con_Graph, int num_nodes = 100, int num_para_per_node = 3, int num_neighbor_radius = 3)
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
