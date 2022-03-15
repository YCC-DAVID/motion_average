#include "MoAve_0.h"




// initialize class specific variables
MoAve_0::MoAve_0()
{



	
}


MoAve_0::~MoAve_0()
{




}
/*******************************************************************************************
direct solver
*******************************************************************************************/
bool MoAve_0::DirectSolver_Run()
{





	return true;
}

/**************************************************************************************
check the parent class for instruction. This make updates for the Node_para_list_buf
covariance for describing the contribution of different parameters of the parameter set.
note, the pair-wise covariance matrix will not be used.
make sure the relax_factor is there.
The P_ij is stored in _Con_Graph
***************************************************************************************/
bool MoAve_0::_LocalSolver(int index_i)
{
	
	// Update through the _Con_Graph, CovMatrix   
	for (int p = index_i; p <= index_i; p ++)
	{
		if (_Ref_Node_Status[p])               // this is the reference point, no updates
			continue;                  

		// aggregate the propagated measurement to this node
		std::vector<double> temvec(_num_para_per_node,0);
		std::vector<double> temWeightVec(_num_para_per_node,0);

		for (int p1 = 0; p1 < (*_Con_Graph)[p].size(); p1 ++)
		{
			/*********************************************************************
			this part can be changed to allow different types of featureless computation
			Node, The _relax_factor is used here in the following, the goal is to update the 
			*********************************************************************/
			
			if (_Ref_Node_Status[(*_Con_Graph)[p][p1].OtherID])  // if this node is linked to reference, give higher weight
			{
				
				for (int p2 = 0; p2 < _num_para_per_node; p2 ++)
				{
					double CurTemWeight = (*_Con_Graph)[p][p1].InvCovMatrix[p2*_num_para_per_node+p2]*(*_Con_Graph)[p][p1].LinkWeight* _Relax_Factor ;
					temvec[p2] += CurTemWeight * (_Node_Para_list[(*_Con_Graph)[p][p1].OtherID][p2] + (*_Con_Graph)[p][p1].RelParaVec[p2]);
				    temWeightVec[p2] += CurTemWeight;
				}


			} // if 
			else
			{
				for (int p2 = 0; p2 < _num_para_per_node; p2 ++)
				{
					double CurTemWeight = (*_Con_Graph)[p][p1].InvCovMatrix[p2*_num_para_per_node+p2]*(*_Con_Graph)[p][p1].LinkWeight;
					temvec[p2] += CurTemWeight * (_Node_Para_list[(*_Con_Graph)[p][p1].OtherID][p2] + (*_Con_Graph)[p][p1].RelParaVec[p2]);
					temWeightVec[p2] += CurTemWeight; 
				}
			} // else

		}// end for p1

		for (int p2 = 0; p2 < _num_para_per_node; p2 ++)
			_Node_Para_list_buf[p][p2] = temvec[p2]/temWeightVec[p2];             // averaging the current updates

	} // end for p (each node) 

	return true;
}

