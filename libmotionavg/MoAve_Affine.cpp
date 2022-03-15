#include "MoAve_Affine.h"




MoAve_Affine::MoAve_Affine()
{





}


MoAve_Affine::~MoAve_Affine()
{





}


bool MoAve_Affine::DirectSolver_Run()
{
	



	return true;
}

/****************************************************************************************************
x2 = Node_para[0] + Node_para[1] * x1 + Node_para[2] * y1;
y2 = Node_para[3] + Node_para[4] * x1 + Node_para[5] * y1;

******************************************************************************************************/


bool MoAve_Affine::_LocalSolver(int index_i)
{

		if (_num_para_per_node != 6)
		{
			std::cout<<"number of parameter per node does not match the problem"<<std::endl;
			return false;
		}

		std::vector<double> temvec(_num_para_per_node,0);
		std::vector<double> temWeightVec(_num_para_per_node,0);

		for (int p1 = 0; p1 < (*_Con_Graph)[index_i].size(); p1 ++)
		{
			/*********************************************************************
			this part can be changed to allow different types of featureless computation
			Node, The _relax_factor is used here in the following, the goal is to update the 
			*********************************************************************/

				double CurTemWeight = 1;
				double Multiplier = 1;

				if (_Ref_Node_Status[(*_Con_Graph)[index_i][p1].OtherID])  //// if this node is linked to reference, give higher weight
					Multiplier = _Relax_Factor;

				// 0
				CurTemWeight = (*_Con_Graph)[index_i][p1].InvCovMatrix[0*_num_para_per_node+0]*(*_Con_Graph)[index_i][p1].LinkWeight* Multiplier;

				temvec[0] +=  _Node_Para_list[index_i][1]*(*_Con_Graph)[index_i][p1].RelParaVec[0]
				+ _Node_Para_list[index_i][2]*(*_Con_Graph)[index_i][p1].RelParaVec[3]
				+ _Node_Para_list[index_i][0];

				temWeightVec[0] += CurTemWeight;

				CurTemWeight = (*_Con_Graph)[index_i][p1].InvCovMatrix[1*_num_para_per_node+1]*(*_Con_Graph)[index_i][p1].LinkWeight* Multiplier ;

				// 1
				temvec[1] +=  _Node_Para_list[index_i][1]*(*_Con_Graph)[index_i][p1].RelParaVec[1]
				+ _Node_Para_list[index_i][2]*(*_Con_Graph)[index_i][p1].RelParaVec[4];

				temWeightVec[1] += CurTemWeight;

				// 2
				CurTemWeight = (*_Con_Graph)[index_i][p1].InvCovMatrix[2*_num_para_per_node+2]*(*_Con_Graph)[index_i][p1].LinkWeight* Multiplier ;

				temvec[2] +=  _Node_Para_list[index_i][1]*(*_Con_Graph)[index_i][p1].RelParaVec[2]
				+ _Node_Para_list[index_i][2]*(*_Con_Graph)[index_i][p1].RelParaVec[5];
				
				temWeightVec[2] += CurTemWeight;


				// 3
				CurTemWeight = (*_Con_Graph)[index_i][p1].InvCovMatrix[3*_num_para_per_node+3]*(*_Con_Graph)[index_i][p1].LinkWeight* Multiplier ;

				temvec[3] +=  _Node_Para_list[index_i][4]*(*_Con_Graph)[index_i][p1].RelParaVec[0]
				+ _Node_Para_list[index_i][5]*(*_Con_Graph)[index_i][p1].RelParaVec[3]
				+ _Node_Para_list[index_i][3];

				temWeightVec[3] += CurTemWeight;

				// 4
				CurTemWeight = (*_Con_Graph)[index_i][p1].InvCovMatrix[4*_num_para_per_node+4]*(*_Con_Graph)[index_i][p1].LinkWeight* Multiplier ;

				temvec[4] +=  _Node_Para_list[index_i][4]*(*_Con_Graph)[index_i][p1].RelParaVec[1]
				+ _Node_Para_list[index_i][5]*(*_Con_Graph)[index_i][p1].RelParaVec[4]
				;
				temWeightVec[4] += CurTemWeight;

				// 5
				CurTemWeight = (*_Con_Graph)[index_i][p1].InvCovMatrix[5*_num_para_per_node+5]*(*_Con_Graph)[index_i][p1].LinkWeight* Multiplier ;

				temvec[5] +=  _Node_Para_list[index_i][4]*(*_Con_Graph)[index_i][p1].RelParaVec[2]
				+ _Node_Para_list[index_i][5]*(*_Con_Graph)[index_i][p1].RelParaVec[5];

				temWeightVec[5] += CurTemWeight;

		}// end for p1

		for (int p2 = 0; p2 < _num_para_per_node; p2 ++)
			_Node_Para_list_buf[index_i][p2] = temvec[p2]/temWeightVec[p2];             // averaging the current updates

	return true;
}