#include <vector>
#include "MoAve_0.h"
#include "MoAve_Affine.h"

int main(int argc, char** argv)
{
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
	M0.Initialize(NodePara, 1, Con_Graph);
	//M0.Set_Reference_Nodes(RefNode);
	M0.IterSolver_Run();

	return true;
}