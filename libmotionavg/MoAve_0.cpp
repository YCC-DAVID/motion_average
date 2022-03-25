#include "MoAve_0.h"

#define USE_CERES_SOLVER

#if defined(USE_CERES_SOLVER)
#include <ceres/problem.h>
#include <ceres/solver.h>
#include <ceres/local_parameterization.h>
#include <ceres/sized_cost_function.h>
#include <ceres/autodiff_cost_function.h>
#include <ceres/loss_function.h>
#endif


// initialize class specific variables
MoAve_0::MoAve_0()
{



	
}


MoAve_0::~MoAve_0()
{




}


#ifdef USE_CERES_SOLVER

struct Translate3WithPriorAD {

	Translate3WithPriorAD(const double* data, const double* weight)
	{
		std::copy_n(data, 3, xfm);
		std::copy_n(weight, 9, P);
	}

	template <typename T>
	bool operator()(const T* v0,
		const T* v1,
		T* residuals) const
	{
		T tem[3];
		for(int i=0;i<3;++i)
			tem[i] = v0[i] - v1[i] - xfm[i];

		for (int i = 0; i < 3; ++i)
		{
			residuals[i] = T(0.);
			for (int j = 0; j < 3; ++j)
				residuals[i] += P[3 * i + j] * tem[j];
		}
		return true;
	}

	static ceres::CostFunction* Create(const double* data, const double* weight)
	{
		return (new ceres::AutoDiffCostFunction<Translate3WithPriorAD, 3, 3, 3>(new Translate3WithPriorAD(data, weight)));
	}
	double xfm[3];
	double P[9]; // 3x3 weight matrix
};

#endif

/*******************************************************************************************
direct solver
*******************************************************************************************/
bool MoAve_0::DirectSolver_Run()
{
	if (_num_para_per_node != 3)
	{
		std::cout << "number of parameter per node does not match the problem" << std::endl;
		return false;
	}
#if defined(USE_CERES_SOLVER)
	// Declare problem
	ceres::Problem problem;

	size_t num_nodes = std::min<size_t>(_Node_Para_list.size(), _Con_Graph->size());


	for (size_t ni = 0; ni < num_nodes; ++ni)
	{
		for (size_t ei = 0; ei < (*_Con_Graph)[ni].size(); ++ei)
		{
			const auto& e = (*_Con_Graph)[ni][ei];

			if (e.RefID != ni)
			{
				std::cout << "Edge.RefID " << e.RefID << " doesn't much Con_Graph block " << ni << std::endl;
				continue;
			}

			ceres::CostFunction* cost_function = Translate3WithPriorAD::Create(e.RelParaVec.data(), e.InvCovMatrix.data());
			//ceres::CostFunction* cost_function = new Translate3WithPriorAD(e.RelParaVec.data(), e.InvCovMatrix.data());
			problem.AddResidualBlock(cost_function,
				nullptr,
				_Node_Para_list[e.RefID].data(),
				_Node_Para_list[e.OtherID].data());
		}
	}
	// Solve problem
	ceres::Solver::Options options;
	options.linear_solver_type = ceres::DENSE_SCHUR;
	options.minimizer_progress_to_stdout = true;
	//options.check_gradients = true;
	options.function_tolerance = _epslon * 1e-4;
	options.gradient_tolerance = _epslon * 1e-4;
	options.parameter_tolerance = _epslon * 1e-4;
	options.max_num_iterations = _num_max_Iter;

	std::vector<int> sspv;
	for (int i = 0; i < _num_para_per_node; ++i)
		sspv.push_back(i);
	auto cssp = new ceres::SubsetParameterization(_num_para_per_node, sspv);
	for (auto rid : _RefNode_List)
	{
		problem.SetParameterization(_Node_Para_list[rid].data(), cssp);
	}

	ceres::Solver::Summary summary;
	ceres::Solve(options, &problem, &summary);
	std::cout << summary.FullReport() << std::endl;
	return true;

#elif defined(USE_ALGLIB)

#endif
	return false;
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

