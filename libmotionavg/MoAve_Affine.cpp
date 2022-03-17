#include "MoAve_Affine.h"

#define USE_CERES_SOLVER

#ifdef USE_CERES_SOLVER
#include <ceres/problem.h>
#include <ceres/solver.h>
#include <ceres/local_parameterization.h>
#include <ceres/autodiff_cost_function.h>
#include <ceres/loss_function.h>
#endif

/*
	* Composition of Affine Transformations
	* lMat = [	l(1)	l(2)	l(0);	rMat = [	r(1)	r(2)	r(0);
	*			l(4)	l(5)	l(3);				r(4)	r(5)	r(3);
	*			0		0		1	]				0		0		1	]
	*
	* cMat = [	c(1)	c(2)	c(0);
	*			c(4)	c(5)	c(3);
	*			0		0		1	]
	*
	*	= [	l(1)*r(1)+l(2)*r(4)	l(1)*r(2)+l(2)*r(5)	l(1)*r(0)+l(2)*r(3)+l(0)	;
	*		l(4)*r(1)+l(5)*r(4)	l(4)*r(2)+l(5)*r(5)	l(4)*r(0)+l(5)*r(3)+l(3)	;
	*		0					0					1							]
	*
	*/
template <typename LT, typename RT, typename CT>
void XfmComposite(const LT* const l, const RT* const r, CT* c)
{
	c[0] = l[1] * r[0] + l[2] * r[3] + l[0];
	c[1] = l[1] * r[1] + l[2] * r[4];
	c[2] = l[1] * r[2] + l[2] * r[5];
	c[3] = l[4] * r[0] + l[5] * r[3] + l[3];
	c[4] = l[4] * r[1] + l[5] * r[4];
	c[5] = l[4] * r[2] + l[5] * r[5];
}

/*
* Inverse of affine transformation
* vMat =
*
* D = (v1*v5 - v2*v4)
*
* inv(vMat) = [	i(1)	i(2)	i(0);
*				i(4)	i(5)	i(3);
*				0		0		1	]
*	=	[  v5/D		-v2/D	-(v0*v5 - v2*v3)/D	;
*			-v4/D	v1/D	(v0*v4 - v1*v3)/D	;
*			0		0		1					]
*/
template <typename IT, typename OT>
void XfmInv(const IT* v, OT* i)
{
	IT D = (v[1] * v[5] - v[2] * v[4]);
	i[0] = -(v[0] * v[5] - v[2] * v[3]) / D;
	i[1] = v[5] / D;
	i[2] = -v[2] / D;
	i[3] = (v[0] * v[4] - v[1] * v[3]) / D;
	i[4] = -v[4] / D;
	i[5] = v[1] / D;
}

template<typename PT, typename IT, typename OT>
void XfmApply(const PT* param, const IT x, const IT y, OT& ox, OT& oy) {
	ox = param[0] + param[1] * x + param[2] * y;
	oy = param[3] + param[4] * x + param[5] * y;
}


MoAve_Affine::MoAve_Affine()
{





}


MoAve_Affine::~MoAve_Affine()
{





}

#ifdef USE_CERES_SOLVER

struct AdditiveWithPrior_PoseError {

	AdditiveWithPrior_PoseError(const double* data, const double* weight)
	{
		std::copy_n(data, 6, xfm);
		std::copy_n(weight, 36, P);
	}
	/*
	* rel_hat = inv(v0)*v1
	* xfm = observed(inv(v0)*reg*v1)
	* r = v0 - v0_hat
	*/
	template <typename T>
	bool operator()(const T* v0,
		const T* v1,
		T* residuals) const
	{
		T inv0[6];
		XfmInv(v0, inv0);
		T rel01[6];
		XfmComposite(inv0, v1, rel01);

		T res[6];
		for (int i = 0; i < 6; ++i)
			res[i] = rel01[i] - xfm[i];

		for (int i = 0; i < 6; ++i)
		{
			residuals[i] = T(0);
			for (int j = 0; j < 6; ++j)
			{
				residuals[i] += res[j] * P[i * 6 + j];
			}
		}
		return true;
	}

	static ceres::CostFunction* Create(const double* data, const double* weight)
	{
		return (new ceres::AutoDiffCostFunction<AdditiveWithPrior_PoseError, 6, 6, 6>(new AdditiveWithPrior_PoseError(data, weight)));
	}
	double xfm[6]; // tx a11 a12 ty a21 a22
	double P[36]; // 6x6 weight matrix
};
#endif

bool MoAve_Affine::DirectSolver_Run()
{
	if (_num_para_per_node != 6)
	{
		std::cout << "number of parameter per node does not match the problem" << std::endl;
		return false;
	}
#ifdef USE_CERES_SOLVER
	// Declare problem
	ceres::Problem problem;

	size_t num_nodes = std::min<size_t>(_Node_Para_list.size(),_Con_Graph->size());


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

			ceres::CostFunction* cost_function = AdditiveWithPrior_PoseError::Create(e.RelParaVec.data(),e.InvCovMatrix.data());
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
	options.function_tolerance = _epslon*1e-4;
	options.gradient_tolerance = _epslon*1e-4;
	options.parameter_tolerance = _epslon*1e-4;
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
#endif

	return false;
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

				const int refID = index_i;
				const int otherID = (*_Con_Graph)[index_i][p1].OtherID;

				std::vector<double>& refNodePara = _Node_Para_list[index_i];
				std::vector<double>& otherNodePara = _Node_Para_list[otherID];
				std::vector<double>& relPara = (*_Con_Graph)[index_i][p1].RelParaVec;
				
				double invRelPara[6], estRef[6];
				XfmInv(relPara.data(), invRelPara);
				XfmComposite( otherNodePara.data(), invRelPara, estRef);
				
				// 0
				CurTemWeight = (*_Con_Graph)[index_i][p1].InvCovMatrix[0*_num_para_per_node+0]*(*_Con_Graph)[index_i][p1].LinkWeight* Multiplier;

				temvec[0] += CurTemWeight * estRef[0];

				temWeightVec[0] += CurTemWeight;

				// 1
				CurTemWeight = (*_Con_Graph)[index_i][p1].InvCovMatrix[1*_num_para_per_node+1]*(*_Con_Graph)[index_i][p1].LinkWeight* Multiplier ;

				temvec[1] += CurTemWeight * estRef[1];

				temWeightVec[1] += CurTemWeight;

				// 2
				CurTemWeight = (*_Con_Graph)[index_i][p1].InvCovMatrix[2*_num_para_per_node+2]*(*_Con_Graph)[index_i][p1].LinkWeight* Multiplier ;

				temvec[2] += CurTemWeight * estRef[2];
				
				temWeightVec[2] += CurTemWeight;


				// 3
				CurTemWeight = (*_Con_Graph)[index_i][p1].InvCovMatrix[3*_num_para_per_node+3]*(*_Con_Graph)[index_i][p1].LinkWeight* Multiplier ;
				
				temvec[3] += CurTemWeight * estRef[3];

				temWeightVec[3] += CurTemWeight;

				// 4
				CurTemWeight = (*_Con_Graph)[index_i][p1].InvCovMatrix[4*_num_para_per_node+4]*(*_Con_Graph)[index_i][p1].LinkWeight* Multiplier ;

				temvec[4] += CurTemWeight * estRef[4];

				temWeightVec[4] += CurTemWeight;

				// 5
				CurTemWeight = (*_Con_Graph)[index_i][p1].InvCovMatrix[5*_num_para_per_node+5]*(*_Con_Graph)[index_i][p1].LinkWeight* Multiplier ;

				temvec[5] += CurTemWeight * estRef[5];

				temWeightVec[5] += CurTemWeight;

		}// end for p1

		for (int p2 = 0; p2 < _num_para_per_node; p2 ++)
			_Node_Para_list_buf[index_i][p2] = temvec[p2]/temWeightVec[p2];             // averaging the current updates

	return true;
}