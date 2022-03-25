#include "MoAve_Affine.h"

#define USE_CERES_SOLVER

#if defined(USE_CERES_SOLVER)
#include <ceres/problem.h>
#include <ceres/solver.h>
#include <ceres/local_parameterization.h>
#include <ceres/sized_cost_function.h>
#include <ceres/autodiff_cost_function.h>
#include <ceres/loss_function.h>
#endif

/*
	* Composition of Affine Transformations
	* lMat = [	l(1)	l(2)	l(0);	rMat = [	v1(1)	v1(2)	v1(0);
	*			l(4)	l(5)	l(3);				v1(4)	v1(5)	v1(3);
	*			0		0		1	]				0		0		1	]
	*
	* cMat = [	c(1)	c(2)	c(0);
	*			c(4)	c(5)	c(3);
	*			0		0		1	]
	*
	*	= [	l(1)*v1(1)+l(2)*v1(4)	l(1)*v1(2)+l(2)*v1(5)	l(1)*v1(0)+l(2)*v1(3)+l(0)	;
	*		l(4)*v1(1)+l(5)*v1(4)	l(4)*v1(2)+l(5)*v1(5)	l(4)*v1(0)+l(5)*v1(3)+l(3)	;
	*		0					0					1							]
	*
	*/
template <typename LT, typename RT, typename CT>
void XfmComposite(const LT* const l, const RT* const v1, CT* c)
{
	c[0] = l[1] * v1[0] + l[2] * v1[3] + l[0];
	c[1] = l[1] * v1[1] + l[2] * v1[4];
	c[2] = l[1] * v1[2] + l[2] * v1[5];
	c[3] = l[4] * v1[0] + l[5] * v1[3] + l[3];
	c[4] = l[4] * v1[1] + l[5] * v1[4];
	c[5] = l[4] * v1[2] + l[5] * v1[5];
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
	* v1 = v0 - v0_hat
	*/
	template <typename T>
	bool operator()(const T* v0,
		const T* v1,
		T* residuals) const
	{
		residuals[0] = P[0] * (-v0[2] * v1[3] / (v0[1] * v0[5] - v0[2] * v0[4]) + v0[5] * v1[0] / (v0[1] * v0[5] - v0[2] * v0[4]) - xfm[0] + (-v0[0] * v0[5] + v0[2] * v0[3]) / (v0[1] * v0[5] - v0[2] * v0[4])) + P[1] * (-v0[2] * v1[4] / (v0[1] * v0[5] - v0[2] * v0[4]) + v0[5] * v1[1] / (v0[1] * v0[5] - v0[2] * v0[4]) - xfm[1]) + P[2] * (-v0[2] * v1[5] / (v0[1] * v0[5] - v0[2] * v0[4]) + v0[5] * v1[2] / (v0[1] * v0[5] - v0[2] * v0[4]) - xfm[2]) + P[3] * (v0[1] * v1[3] / (v0[1] * v0[5] - v0[2] * v0[4]) - v0[4] * v1[0] / (v0[1] * v0[5] - v0[2] * v0[4]) - xfm[3] + (v0[0] * v0[4] - v0[1] * v0[3]) / (v0[1] * v0[5] - v0[2] * v0[4])) + P[4] * (v0[1] * v1[4] / (v0[1] * v0[5] - v0[2] * v0[4]) - v0[4] * v1[1] / (v0[1] * v0[5] - v0[2] * v0[4]) - xfm[4]) + P[5] * (v0[1] * v1[5] / (v0[1] * v0[5] - v0[2] * v0[4]) - v0[4] * v1[2] / (v0[1] * v0[5] - v0[2] * v0[4]) - xfm[5]);
		residuals[1] = P[10] * (v0[1] * v1[4] / (v0[1] * v0[5] - v0[2] * v0[4]) - v0[4] * v1[1] / (v0[1] * v0[5] - v0[2] * v0[4]) - xfm[4]) + P[11] * (v0[1] * v1[5] / (v0[1] * v0[5] - v0[2] * v0[4]) - v0[4] * v1[2] / (v0[1] * v0[5] - v0[2] * v0[4]) - xfm[5]) + P[6] * (-v0[2] * v1[3] / (v0[1] * v0[5] - v0[2] * v0[4]) + v0[5] * v1[0] / (v0[1] * v0[5] - v0[2] * v0[4]) - xfm[0] + (-v0[0] * v0[5] + v0[2] * v0[3]) / (v0[1] * v0[5] - v0[2] * v0[4])) + P[7] * (-v0[2] * v1[4] / (v0[1] * v0[5] - v0[2] * v0[4]) + v0[5] * v1[1] / (v0[1] * v0[5] - v0[2] * v0[4]) - xfm[1]) + P[8] * (-v0[2] * v1[5] / (v0[1] * v0[5] - v0[2] * v0[4]) + v0[5] * v1[2] / (v0[1] * v0[5] - v0[2] * v0[4]) - xfm[2]) + P[9] * (v0[1] * v1[3] / (v0[1] * v0[5] - v0[2] * v0[4]) - v0[4] * v1[0] / (v0[1] * v0[5] - v0[2] * v0[4]) - xfm[3] + (v0[0] * v0[4] - v0[1] * v0[3]) / (v0[1] * v0[5] - v0[2] * v0[4]));
		residuals[2] = P[12] * (-v0[2] * v1[3] / (v0[1] * v0[5] - v0[2] * v0[4]) + v0[5] * v1[0] / (v0[1] * v0[5] - v0[2] * v0[4]) - xfm[0] + (-v0[0] * v0[5] + v0[2] * v0[3]) / (v0[1] * v0[5] - v0[2] * v0[4])) + P[13] * (-v0[2] * v1[4] / (v0[1] * v0[5] - v0[2] * v0[4]) + v0[5] * v1[1] / (v0[1] * v0[5] - v0[2] * v0[4]) - xfm[1]) + P[14] * (-v0[2] * v1[5] / (v0[1] * v0[5] - v0[2] * v0[4]) + v0[5] * v1[2] / (v0[1] * v0[5] - v0[2] * v0[4]) - xfm[2]) + P[15] * (v0[1] * v1[3] / (v0[1] * v0[5] - v0[2] * v0[4]) - v0[4] * v1[0] / (v0[1] * v0[5] - v0[2] * v0[4]) - xfm[3] + (v0[0] * v0[4] - v0[1] * v0[3]) / (v0[1] * v0[5] - v0[2] * v0[4])) + P[16] * (v0[1] * v1[4] / (v0[1] * v0[5] - v0[2] * v0[4]) - v0[4] * v1[1] / (v0[1] * v0[5] - v0[2] * v0[4]) - xfm[4]) + P[17] * (v0[1] * v1[5] / (v0[1] * v0[5] - v0[2] * v0[4]) - v0[4] * v1[2] / (v0[1] * v0[5] - v0[2] * v0[4]) - xfm[5]);
		residuals[3] = P[18] * (-v0[2] * v1[3] / (v0[1] * v0[5] - v0[2] * v0[4]) + v0[5] * v1[0] / (v0[1] * v0[5] - v0[2] * v0[4]) - xfm[0] + (-v0[0] * v0[5] + v0[2] * v0[3]) / (v0[1] * v0[5] - v0[2] * v0[4])) + P[19] * (-v0[2] * v1[4] / (v0[1] * v0[5] - v0[2] * v0[4]) + v0[5] * v1[1] / (v0[1] * v0[5] - v0[2] * v0[4]) - xfm[1]) + P[20] * (-v0[2] * v1[5] / (v0[1] * v0[5] - v0[2] * v0[4]) + v0[5] * v1[2] / (v0[1] * v0[5] - v0[2] * v0[4]) - xfm[2]) + P[21] * (v0[1] * v1[3] / (v0[1] * v0[5] - v0[2] * v0[4]) - v0[4] * v1[0] / (v0[1] * v0[5] - v0[2] * v0[4]) - xfm[3] + (v0[0] * v0[4] - v0[1] * v0[3]) / (v0[1] * v0[5] - v0[2] * v0[4])) + P[22] * (v0[1] * v1[4] / (v0[1] * v0[5] - v0[2] * v0[4]) - v0[4] * v1[1] / (v0[1] * v0[5] - v0[2] * v0[4]) - xfm[4]) + P[23] * (v0[1] * v1[5] / (v0[1] * v0[5] - v0[2] * v0[4]) - v0[4] * v1[2] / (v0[1] * v0[5] - v0[2] * v0[4]) - xfm[5]);
		residuals[4] = P[24] * (-v0[2] * v1[3] / (v0[1] * v0[5] - v0[2] * v0[4]) + v0[5] * v1[0] / (v0[1] * v0[5] - v0[2] * v0[4]) - xfm[0] + (-v0[0] * v0[5] + v0[2] * v0[3]) / (v0[1] * v0[5] - v0[2] * v0[4])) + P[25] * (-v0[2] * v1[4] / (v0[1] * v0[5] - v0[2] * v0[4]) + v0[5] * v1[1] / (v0[1] * v0[5] - v0[2] * v0[4]) - xfm[1]) + P[26] * (-v0[2] * v1[5] / (v0[1] * v0[5] - v0[2] * v0[4]) + v0[5] * v1[2] / (v0[1] * v0[5] - v0[2] * v0[4]) - xfm[2]) + P[27] * (v0[1] * v1[3] / (v0[1] * v0[5] - v0[2] * v0[4]) - v0[4] * v1[0] / (v0[1] * v0[5] - v0[2] * v0[4]) - xfm[3] + (v0[0] * v0[4] - v0[1] * v0[3]) / (v0[1] * v0[5] - v0[2] * v0[4])) + P[28] * (v0[1] * v1[4] / (v0[1] * v0[5] - v0[2] * v0[4]) - v0[4] * v1[1] / (v0[1] * v0[5] - v0[2] * v0[4]) - xfm[4]) + P[29] * (v0[1] * v1[5] / (v0[1] * v0[5] - v0[2] * v0[4]) - v0[4] * v1[2] / (v0[1] * v0[5] - v0[2] * v0[4]) - xfm[5]);
		residuals[5] = P[30] * (-v0[2] * v1[3] / (v0[1] * v0[5] - v0[2] * v0[4]) + v0[5] * v1[0] / (v0[1] * v0[5] - v0[2] * v0[4]) - xfm[0] + (-v0[0] * v0[5] + v0[2] * v0[3]) / (v0[1] * v0[5] - v0[2] * v0[4])) + P[31] * (-v0[2] * v1[4] / (v0[1] * v0[5] - v0[2] * v0[4]) + v0[5] * v1[1] / (v0[1] * v0[5] - v0[2] * v0[4]) - xfm[1]) + P[32] * (-v0[2] * v1[5] / (v0[1] * v0[5] - v0[2] * v0[4]) + v0[5] * v1[2] / (v0[1] * v0[5] - v0[2] * v0[4]) - xfm[2]) + P[33] * (v0[1] * v1[3] / (v0[1] * v0[5] - v0[2] * v0[4]) - v0[4] * v1[0] / (v0[1] * v0[5] - v0[2] * v0[4]) - xfm[3] + (v0[0] * v0[4] - v0[1] * v0[3]) / (v0[1] * v0[5] - v0[2] * v0[4])) + P[34] * (v0[1] * v1[4] / (v0[1] * v0[5] - v0[2] * v0[4]) - v0[4] * v1[1] / (v0[1] * v0[5] - v0[2] * v0[4]) - xfm[4]) + P[35] * (v0[1] * v1[5] / (v0[1] * v0[5] - v0[2] * v0[4]) - v0[4] * v1[2] / (v0[1] * v0[5] - v0[2] * v0[4]) - xfm[5]);

		return true;
	}

	static ceres::CostFunction* Create(const double* data, const double* weight)
	{
		return (new ceres::AutoDiffCostFunction<AdditiveWithPrior_PoseError, 6, 6, 6>(new AdditiveWithPrior_PoseError(data, weight)));
	}
	double xfm[6]; // tx a11 a12 ty a21 a22
	double P[36]; // 6x6 weight matrix
};

class Affine2WithPrior_Additive_Analytic : public ceres::SizedCostFunction<6, 6, 6> {
public:
	Affine2WithPrior_Additive_Analytic(const double* data, const double* weight)
	{
		std::copy_n(data, 6, xfm);
		std::copy_n(weight, 36, P);
	}
	virtual ~Affine2WithPrior_Additive_Analytic() {}
	virtual bool Evaluate(double const* const* parameters,
		double* residuals,
		double** jacobians) const {

		double const* v0 = parameters[0];
		double const* v1 = parameters[1];

		residuals[0] = (-P[0] * (v0[0] * v0[5] - v0[2] * v0[3] + v0[2] * v1[3] - v0[5] * v1[0] + xfm[0] * (v0[1] * v0[5] - v0[2] * v0[4])) - P[1] * (v0[2] * v1[4] - v0[5] * v1[1] + xfm[1] * (v0[1] * v0[5] - v0[2] * v0[4])) - P[2] * (v0[2] * v1[5] - v0[5] * v1[2] + xfm[2] * (v0[1] * v0[5] - v0[2] * v0[4])) + P[3] * (v0[0] * v0[4] - v0[1] * v0[3] + v0[1] * v1[3] - v0[4] * v1[0] - xfm[3] * (v0[1] * v0[5] - v0[2] * v0[4])) - P[4] * (-v0[1] * v1[4] + v0[4] * v1[1] + xfm[4] * (v0[1] * v0[5] - v0[2] * v0[4])) - P[5] * (-v0[1] * v1[5] + v0[4] * v1[2] + xfm[5] * (v0[1] * v0[5] - v0[2] * v0[4]))) / (v0[1] * v0[5] - v0[2] * v0[4]);
		residuals[1] = (-P[10] * (-v0[1] * v1[4] + v0[4] * v1[1] + xfm[4] * (v0[1] * v0[5] - v0[2] * v0[4])) - P[11] * (-v0[1] * v1[5] + v0[4] * v1[2] + xfm[5] * (v0[1] * v0[5] - v0[2] * v0[4])) - P[6] * (v0[0] * v0[5] - v0[2] * v0[3] + v0[2] * v1[3] - v0[5] * v1[0] + xfm[0] * (v0[1] * v0[5] - v0[2] * v0[4])) - P[7] * (v0[2] * v1[4] - v0[5] * v1[1] + xfm[1] * (v0[1] * v0[5] - v0[2] * v0[4])) - P[8] * (v0[2] * v1[5] - v0[5] * v1[2] + xfm[2] * (v0[1] * v0[5] - v0[2] * v0[4])) + P[9] * (v0[0] * v0[4] - v0[1] * v0[3] + v0[1] * v1[3] - v0[4] * v1[0] - xfm[3] * (v0[1] * v0[5] - v0[2] * v0[4]))) / (v0[1] * v0[5] - v0[2] * v0[4]);
		residuals[2] = (-P[12] * (v0[0] * v0[5] - v0[2] * v0[3] + v0[2] * v1[3] - v0[5] * v1[0] + xfm[0] * (v0[1] * v0[5] - v0[2] * v0[4])) - P[13] * (v0[2] * v1[4] - v0[5] * v1[1] + xfm[1] * (v0[1] * v0[5] - v0[2] * v0[4])) - P[14] * (v0[2] * v1[5] - v0[5] * v1[2] + xfm[2] * (v0[1] * v0[5] - v0[2] * v0[4])) + P[15] * (v0[0] * v0[4] - v0[1] * v0[3] + v0[1] * v1[3] - v0[4] * v1[0] - xfm[3] * (v0[1] * v0[5] - v0[2] * v0[4])) - P[16] * (-v0[1] * v1[4] + v0[4] * v1[1] + xfm[4] * (v0[1] * v0[5] - v0[2] * v0[4])) - P[17] * (-v0[1] * v1[5] + v0[4] * v1[2] + xfm[5] * (v0[1] * v0[5] - v0[2] * v0[4]))) / (v0[1] * v0[5] - v0[2] * v0[4]);
		residuals[3] = (-P[18] * (v0[0] * v0[5] - v0[2] * v0[3] + v0[2] * v1[3] - v0[5] * v1[0] + xfm[0] * (v0[1] * v0[5] - v0[2] * v0[4])) - P[19] * (v0[2] * v1[4] - v0[5] * v1[1] + xfm[1] * (v0[1] * v0[5] - v0[2] * v0[4])) - P[20] * (v0[2] * v1[5] - v0[5] * v1[2] + xfm[2] * (v0[1] * v0[5] - v0[2] * v0[4])) + P[21] * (v0[0] * v0[4] - v0[1] * v0[3] + v0[1] * v1[3] - v0[4] * v1[0] - xfm[3] * (v0[1] * v0[5] - v0[2] * v0[4])) - P[22] * (-v0[1] * v1[4] + v0[4] * v1[1] + xfm[4] * (v0[1] * v0[5] - v0[2] * v0[4])) - P[23] * (-v0[1] * v1[5] + v0[4] * v1[2] + xfm[5] * (v0[1] * v0[5] - v0[2] * v0[4]))) / (v0[1] * v0[5] - v0[2] * v0[4]);
		residuals[4] = (-P[24] * (v0[0] * v0[5] - v0[2] * v0[3] + v0[2] * v1[3] - v0[5] * v1[0] + xfm[0] * (v0[1] * v0[5] - v0[2] * v0[4])) - P[25] * (v0[2] * v1[4] - v0[5] * v1[1] + xfm[1] * (v0[1] * v0[5] - v0[2] * v0[4])) - P[26] * (v0[2] * v1[5] - v0[5] * v1[2] + xfm[2] * (v0[1] * v0[5] - v0[2] * v0[4])) + P[27] * (v0[0] * v0[4] - v0[1] * v0[3] + v0[1] * v1[3] - v0[4] * v1[0] - xfm[3] * (v0[1] * v0[5] - v0[2] * v0[4])) - P[28] * (-v0[1] * v1[4] + v0[4] * v1[1] + xfm[4] * (v0[1] * v0[5] - v0[2] * v0[4])) - P[29] * (-v0[1] * v1[5] + v0[4] * v1[2] + xfm[5] * (v0[1] * v0[5] - v0[2] * v0[4]))) / (v0[1] * v0[5] - v0[2] * v0[4]);
		residuals[5] = (-P[30] * (v0[0] * v0[5] - v0[2] * v0[3] + v0[2] * v1[3] - v0[5] * v1[0] + xfm[0] * (v0[1] * v0[5] - v0[2] * v0[4])) - P[31] * (v0[2] * v1[4] - v0[5] * v1[1] + xfm[1] * (v0[1] * v0[5] - v0[2] * v0[4])) - P[32] * (v0[2] * v1[5] - v0[5] * v1[2] + xfm[2] * (v0[1] * v0[5] - v0[2] * v0[4])) + P[33] * (v0[0] * v0[4] - v0[1] * v0[3] + v0[1] * v1[3] - v0[4] * v1[0] - xfm[3] * (v0[1] * v0[5] - v0[2] * v0[4])) - P[34] * (-v0[1] * v1[4] + v0[4] * v1[1] + xfm[4] * (v0[1] * v0[5] - v0[2] * v0[4])) - P[35] * (-v0[1] * v1[5] + v0[4] * v1[2] + xfm[5] * (v0[1] * v0[5] - v0[2] * v0[4]))) / (v0[1] * v0[5] - v0[2] * v0[4]);

		if (!jacobians) return true;
		
		// Jacobian
		if (jacobians[0]) {
			jacobians[0][0] = (-P[0] * v0[5] + P[3] * v0[4]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[0][1] = (-P[3] * (v0[5] * (v0[0] * v0[4] - v0[1] * v0[3] + v0[1] * v1[3] - v0[4] * v1[0]) + (v0[3] - v1[3]) * (v0[1] * v0[5] - v0[2] * v0[4])) + P[4] * (v0[5] * (-v0[1] * v1[4] + v0[4] * v1[1]) + v1[4] * (v0[1] * v0[5] - v0[2] * v0[4])) + P[5] * (v0[5] * (-v0[1] * v1[5] + v0[4] * v1[2]) + v1[5] * (v0[1] * v0[5] - v0[2] * v0[4])) + v0[5] * (P[0] * (v0[0] * v0[5] - v0[2] * v0[3] + v0[2] * v1[3] - v0[5] * v1[0]) + P[1] * (v0[2] * v1[4] - v0[5] * v1[1]) + P[2] * (v0[2] * v1[5] - v0[5] * v1[2]))) / std::pow(v0[1] * v0[5] - v0[2] * v0[4], 2);
			jacobians[0][2] = (-P[0] * (v0[4] * (v0[0] * v0[5] - v0[2] * v0[3] + v0[2] * v1[3] - v0[5] * v1[0]) + (-v0[3] + v1[3]) * (v0[1] * v0[5] - v0[2] * v0[4])) - P[1] * (v0[4] * (v0[2] * v1[4] - v0[5] * v1[1]) + v1[4] * (v0[1] * v0[5] - v0[2] * v0[4])) - P[2] * (v0[4] * (v0[2] * v1[5] - v0[5] * v1[2]) + v1[5] * (v0[1] * v0[5] - v0[2] * v0[4])) + v0[4] * (P[3] * (v0[0] * v0[4] - v0[1] * v0[3] + v0[1] * v1[3] - v0[4] * v1[0]) + P[4] * (v0[1] * v1[4] - v0[4] * v1[1]) + P[5] * (v0[1] * v1[5] - v0[4] * v1[2]))) / std::pow(v0[1] * v0[5] - v0[2] * v0[4], 2);
			jacobians[0][3] = (P[0] * v0[2] - P[3] * v0[1]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[0][4] = (P[3] * (v0[2] * (v0[0] * v0[4] - v0[1] * v0[3] + v0[1] * v1[3] - v0[4] * v1[0]) + (v0[0] - v1[0]) * (v0[1] * v0[5] - v0[2] * v0[4])) - P[4] * (v0[2] * (-v0[1] * v1[4] + v0[4] * v1[1]) + v1[1] * (v0[1] * v0[5] - v0[2] * v0[4])) - P[5] * (v0[2] * (-v0[1] * v1[5] + v0[4] * v1[2]) + v1[2] * (v0[1] * v0[5] - v0[2] * v0[4])) - v0[2] * (P[0] * (v0[0] * v0[5] - v0[2] * v0[3] + v0[2] * v1[3] - v0[5] * v1[0]) + P[1] * (v0[2] * v1[4] - v0[5] * v1[1]) + P[2] * (v0[2] * v1[5] - v0[5] * v1[2]))) / std::pow(v0[1] * v0[5] - v0[2] * v0[4], 2);
			jacobians[0][5] = (P[0] * (v0[1] * (v0[0] * v0[5] - v0[2] * v0[3] + v0[2] * v1[3] - v0[5] * v1[0]) + (-v0[0] + v1[0]) * (v0[1] * v0[5] - v0[2] * v0[4])) + P[1] * (v0[1] * (v0[2] * v1[4] - v0[5] * v1[1]) + v1[1] * (v0[1] * v0[5] - v0[2] * v0[4])) + P[2] * (v0[1] * (v0[2] * v1[5] - v0[5] * v1[2]) + v1[2] * (v0[1] * v0[5] - v0[2] * v0[4])) - v0[1] * (P[3] * (v0[0] * v0[4] - v0[1] * v0[3] + v0[1] * v1[3] - v0[4] * v1[0]) + P[4] * (v0[1] * v1[4] - v0[4] * v1[1]) + P[5] * (v0[1] * v1[5] - v0[4] * v1[2]))) / std::pow(v0[1] * v0[5] - v0[2] * v0[4], 2);
			jacobians[0][6] = (-P[6] * v0[5] + P[9] * v0[4]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[0][7] = (P[10] * (v0[5] * (-v0[1] * v1[4] + v0[4] * v1[1]) + v1[4] * (v0[1] * v0[5] - v0[2] * v0[4])) + P[11] * (v0[5] * (-v0[1] * v1[5] + v0[4] * v1[2]) + v1[5] * (v0[1] * v0[5] - v0[2] * v0[4])) - P[9] * (v0[5] * (v0[0] * v0[4] - v0[1] * v0[3] + v0[1] * v1[3] - v0[4] * v1[0]) + (v0[3] - v1[3]) * (v0[1] * v0[5] - v0[2] * v0[4])) + v0[5] * (P[6] * (v0[0] * v0[5] - v0[2] * v0[3] + v0[2] * v1[3] - v0[5] * v1[0]) + P[7] * (v0[2] * v1[4] - v0[5] * v1[1]) + P[8] * (v0[2] * v1[5] - v0[5] * v1[2]))) / std::pow(v0[1] * v0[5] - v0[2] * v0[4], 2);
			jacobians[0][8] = (-P[6] * (v0[4] * (v0[0] * v0[5] - v0[2] * v0[3] + v0[2] * v1[3] - v0[5] * v1[0]) + (-v0[3] + v1[3]) * (v0[1] * v0[5] - v0[2] * v0[4])) - P[7] * (v0[4] * (v0[2] * v1[4] - v0[5] * v1[1]) + v1[4] * (v0[1] * v0[5] - v0[2] * v0[4])) - P[8] * (v0[4] * (v0[2] * v1[5] - v0[5] * v1[2]) + v1[5] * (v0[1] * v0[5] - v0[2] * v0[4])) + v0[4] * (P[10] * (v0[1] * v1[4] - v0[4] * v1[1]) + P[11] * (v0[1] * v1[5] - v0[4] * v1[2]) + P[9] * (v0[0] * v0[4] - v0[1] * v0[3] + v0[1] * v1[3] - v0[4] * v1[0]))) / std::pow(v0[1] * v0[5] - v0[2] * v0[4], 2);
			jacobians[0][9] = (P[6] * v0[2] - P[9] * v0[1]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[0][10] = (-P[10] * (v0[2] * (-v0[1] * v1[4] + v0[4] * v1[1]) + v1[1] * (v0[1] * v0[5] - v0[2] * v0[4])) - P[11] * (v0[2] * (-v0[1] * v1[5] + v0[4] * v1[2]) + v1[2] * (v0[1] * v0[5] - v0[2] * v0[4])) + P[9] * (v0[2] * (v0[0] * v0[4] - v0[1] * v0[3] + v0[1] * v1[3] - v0[4] * v1[0]) + (v0[0] - v1[0]) * (v0[1] * v0[5] - v0[2] * v0[4])) - v0[2] * (P[6] * (v0[0] * v0[5] - v0[2] * v0[3] + v0[2] * v1[3] - v0[5] * v1[0]) + P[7] * (v0[2] * v1[4] - v0[5] * v1[1]) + P[8] * (v0[2] * v1[5] - v0[5] * v1[2]))) / std::pow(v0[1] * v0[5] - v0[2] * v0[4], 2);
			jacobians[0][11] = (P[6] * (v0[1] * (v0[0] * v0[5] - v0[2] * v0[3] + v0[2] * v1[3] - v0[5] * v1[0]) + (-v0[0] + v1[0]) * (v0[1] * v0[5] - v0[2] * v0[4])) + P[7] * (v0[1] * (v0[2] * v1[4] - v0[5] * v1[1]) + v1[1] * (v0[1] * v0[5] - v0[2] * v0[4])) + P[8] * (v0[1] * (v0[2] * v1[5] - v0[5] * v1[2]) + v1[2] * (v0[1] * v0[5] - v0[2] * v0[4])) - v0[1] * (P[10] * (v0[1] * v1[4] - v0[4] * v1[1]) + P[11] * (v0[1] * v1[5] - v0[4] * v1[2]) + P[9] * (v0[0] * v0[4] - v0[1] * v0[3] + v0[1] * v1[3] - v0[4] * v1[0]))) / std::pow(v0[1] * v0[5] - v0[2] * v0[4], 2);
			jacobians[0][12] = (-P[12] * v0[5] + P[15] * v0[4]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[0][13] = (-P[15] * (v0[5] * (v0[0] * v0[4] - v0[1] * v0[3] + v0[1] * v1[3] - v0[4] * v1[0]) + (v0[3] - v1[3]) * (v0[1] * v0[5] - v0[2] * v0[4])) + P[16] * (v0[5] * (-v0[1] * v1[4] + v0[4] * v1[1]) + v1[4] * (v0[1] * v0[5] - v0[2] * v0[4])) + P[17] * (v0[5] * (-v0[1] * v1[5] + v0[4] * v1[2]) + v1[5] * (v0[1] * v0[5] - v0[2] * v0[4])) + v0[5] * (P[12] * (v0[0] * v0[5] - v0[2] * v0[3] + v0[2] * v1[3] - v0[5] * v1[0]) + P[13] * (v0[2] * v1[4] - v0[5] * v1[1]) + P[14] * (v0[2] * v1[5] - v0[5] * v1[2]))) / std::pow(v0[1] * v0[5] - v0[2] * v0[4], 2);
			jacobians[0][14] = (-P[12] * (v0[4] * (v0[0] * v0[5] - v0[2] * v0[3] + v0[2] * v1[3] - v0[5] * v1[0]) + (-v0[3] + v1[3]) * (v0[1] * v0[5] - v0[2] * v0[4])) - P[13] * (v0[4] * (v0[2] * v1[4] - v0[5] * v1[1]) + v1[4] * (v0[1] * v0[5] - v0[2] * v0[4])) - P[14] * (v0[4] * (v0[2] * v1[5] - v0[5] * v1[2]) + v1[5] * (v0[1] * v0[5] - v0[2] * v0[4])) + v0[4] * (P[15] * (v0[0] * v0[4] - v0[1] * v0[3] + v0[1] * v1[3] - v0[4] * v1[0]) + P[16] * (v0[1] * v1[4] - v0[4] * v1[1]) + P[17] * (v0[1] * v1[5] - v0[4] * v1[2]))) / std::pow(v0[1] * v0[5] - v0[2] * v0[4], 2);
			jacobians[0][15] = (P[12] * v0[2] - P[15] * v0[1]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[0][16] = (P[15] * (v0[2] * (v0[0] * v0[4] - v0[1] * v0[3] + v0[1] * v1[3] - v0[4] * v1[0]) + (v0[0] - v1[0]) * (v0[1] * v0[5] - v0[2] * v0[4])) - P[16] * (v0[2] * (-v0[1] * v1[4] + v0[4] * v1[1]) + v1[1] * (v0[1] * v0[5] - v0[2] * v0[4])) - P[17] * (v0[2] * (-v0[1] * v1[5] + v0[4] * v1[2]) + v1[2] * (v0[1] * v0[5] - v0[2] * v0[4])) - v0[2] * (P[12] * (v0[0] * v0[5] - v0[2] * v0[3] + v0[2] * v1[3] - v0[5] * v1[0]) + P[13] * (v0[2] * v1[4] - v0[5] * v1[1]) + P[14] * (v0[2] * v1[5] - v0[5] * v1[2]))) / std::pow(v0[1] * v0[5] - v0[2] * v0[4], 2);
			jacobians[0][17] = (P[12] * (v0[1] * (v0[0] * v0[5] - v0[2] * v0[3] + v0[2] * v1[3] - v0[5] * v1[0]) + (-v0[0] + v1[0]) * (v0[1] * v0[5] - v0[2] * v0[4])) + P[13] * (v0[1] * (v0[2] * v1[4] - v0[5] * v1[1]) + v1[1] * (v0[1] * v0[5] - v0[2] * v0[4])) + P[14] * (v0[1] * (v0[2] * v1[5] - v0[5] * v1[2]) + v1[2] * (v0[1] * v0[5] - v0[2] * v0[4])) - v0[1] * (P[15] * (v0[0] * v0[4] - v0[1] * v0[3] + v0[1] * v1[3] - v0[4] * v1[0]) + P[16] * (v0[1] * v1[4] - v0[4] * v1[1]) + P[17] * (v0[1] * v1[5] - v0[4] * v1[2]))) / std::pow(v0[1] * v0[5] - v0[2] * v0[4], 2);
			jacobians[0][18] = (-P[18] * v0[5] + P[21] * v0[4]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[0][19] = (-P[21] * (v0[5] * (v0[0] * v0[4] - v0[1] * v0[3] + v0[1] * v1[3] - v0[4] * v1[0]) + (v0[3] - v1[3]) * (v0[1] * v0[5] - v0[2] * v0[4])) + P[22] * (v0[5] * (-v0[1] * v1[4] + v0[4] * v1[1]) + v1[4] * (v0[1] * v0[5] - v0[2] * v0[4])) + P[23] * (v0[5] * (-v0[1] * v1[5] + v0[4] * v1[2]) + v1[5] * (v0[1] * v0[5] - v0[2] * v0[4])) + v0[5] * (P[18] * (v0[0] * v0[5] - v0[2] * v0[3] + v0[2] * v1[3] - v0[5] * v1[0]) + P[19] * (v0[2] * v1[4] - v0[5] * v1[1]) + P[20] * (v0[2] * v1[5] - v0[5] * v1[2]))) / std::pow(v0[1] * v0[5] - v0[2] * v0[4], 2);
			jacobians[0][20] = (-P[18] * (v0[4] * (v0[0] * v0[5] - v0[2] * v0[3] + v0[2] * v1[3] - v0[5] * v1[0]) + (-v0[3] + v1[3]) * (v0[1] * v0[5] - v0[2] * v0[4])) - P[19] * (v0[4] * (v0[2] * v1[4] - v0[5] * v1[1]) + v1[4] * (v0[1] * v0[5] - v0[2] * v0[4])) - P[20] * (v0[4] * (v0[2] * v1[5] - v0[5] * v1[2]) + v1[5] * (v0[1] * v0[5] - v0[2] * v0[4])) + v0[4] * (P[21] * (v0[0] * v0[4] - v0[1] * v0[3] + v0[1] * v1[3] - v0[4] * v1[0]) + P[22] * (v0[1] * v1[4] - v0[4] * v1[1]) + P[23] * (v0[1] * v1[5] - v0[4] * v1[2]))) / std::pow(v0[1] * v0[5] - v0[2] * v0[4], 2);
			jacobians[0][21] = (P[18] * v0[2] - P[21] * v0[1]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[0][22] = (P[21] * (v0[2] * (v0[0] * v0[4] - v0[1] * v0[3] + v0[1] * v1[3] - v0[4] * v1[0]) + (v0[0] - v1[0]) * (v0[1] * v0[5] - v0[2] * v0[4])) - P[22] * (v0[2] * (-v0[1] * v1[4] + v0[4] * v1[1]) + v1[1] * (v0[1] * v0[5] - v0[2] * v0[4])) - P[23] * (v0[2] * (-v0[1] * v1[5] + v0[4] * v1[2]) + v1[2] * (v0[1] * v0[5] - v0[2] * v0[4])) - v0[2] * (P[18] * (v0[0] * v0[5] - v0[2] * v0[3] + v0[2] * v1[3] - v0[5] * v1[0]) + P[19] * (v0[2] * v1[4] - v0[5] * v1[1]) + P[20] * (v0[2] * v1[5] - v0[5] * v1[2]))) / std::pow(v0[1] * v0[5] - v0[2] * v0[4], 2);
			jacobians[0][23] = (P[18] * (v0[1] * (v0[0] * v0[5] - v0[2] * v0[3] + v0[2] * v1[3] - v0[5] * v1[0]) + (-v0[0] + v1[0]) * (v0[1] * v0[5] - v0[2] * v0[4])) + P[19] * (v0[1] * (v0[2] * v1[4] - v0[5] * v1[1]) + v1[1] * (v0[1] * v0[5] - v0[2] * v0[4])) + P[20] * (v0[1] * (v0[2] * v1[5] - v0[5] * v1[2]) + v1[2] * (v0[1] * v0[5] - v0[2] * v0[4])) - v0[1] * (P[21] * (v0[0] * v0[4] - v0[1] * v0[3] + v0[1] * v1[3] - v0[4] * v1[0]) + P[22] * (v0[1] * v1[4] - v0[4] * v1[1]) + P[23] * (v0[1] * v1[5] - v0[4] * v1[2]))) / std::pow(v0[1] * v0[5] - v0[2] * v0[4], 2);
			jacobians[0][24] = (-P[24] * v0[5] + P[27] * v0[4]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[0][25] = (-P[27] * (v0[5] * (v0[0] * v0[4] - v0[1] * v0[3] + v0[1] * v1[3] - v0[4] * v1[0]) + (v0[3] - v1[3]) * (v0[1] * v0[5] - v0[2] * v0[4])) + P[28] * (v0[5] * (-v0[1] * v1[4] + v0[4] * v1[1]) + v1[4] * (v0[1] * v0[5] - v0[2] * v0[4])) + P[29] * (v0[5] * (-v0[1] * v1[5] + v0[4] * v1[2]) + v1[5] * (v0[1] * v0[5] - v0[2] * v0[4])) + v0[5] * (P[24] * (v0[0] * v0[5] - v0[2] * v0[3] + v0[2] * v1[3] - v0[5] * v1[0]) + P[25] * (v0[2] * v1[4] - v0[5] * v1[1]) + P[26] * (v0[2] * v1[5] - v0[5] * v1[2]))) / std::pow(v0[1] * v0[5] - v0[2] * v0[4], 2);
			jacobians[0][26] = (-P[24] * (v0[4] * (v0[0] * v0[5] - v0[2] * v0[3] + v0[2] * v1[3] - v0[5] * v1[0]) + (-v0[3] + v1[3]) * (v0[1] * v0[5] - v0[2] * v0[4])) - P[25] * (v0[4] * (v0[2] * v1[4] - v0[5] * v1[1]) + v1[4] * (v0[1] * v0[5] - v0[2] * v0[4])) - P[26] * (v0[4] * (v0[2] * v1[5] - v0[5] * v1[2]) + v1[5] * (v0[1] * v0[5] - v0[2] * v0[4])) + v0[4] * (P[27] * (v0[0] * v0[4] - v0[1] * v0[3] + v0[1] * v1[3] - v0[4] * v1[0]) + P[28] * (v0[1] * v1[4] - v0[4] * v1[1]) + P[29] * (v0[1] * v1[5] - v0[4] * v1[2]))) / std::pow(v0[1] * v0[5] - v0[2] * v0[4], 2);
			jacobians[0][27] = (P[24] * v0[2] - P[27] * v0[1]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[0][28] = (P[27] * (v0[2] * (v0[0] * v0[4] - v0[1] * v0[3] + v0[1] * v1[3] - v0[4] * v1[0]) + (v0[0] - v1[0]) * (v0[1] * v0[5] - v0[2] * v0[4])) - P[28] * (v0[2] * (-v0[1] * v1[4] + v0[4] * v1[1]) + v1[1] * (v0[1] * v0[5] - v0[2] * v0[4])) - P[29] * (v0[2] * (-v0[1] * v1[5] + v0[4] * v1[2]) + v1[2] * (v0[1] * v0[5] - v0[2] * v0[4])) - v0[2] * (P[24] * (v0[0] * v0[5] - v0[2] * v0[3] + v0[2] * v1[3] - v0[5] * v1[0]) + P[25] * (v0[2] * v1[4] - v0[5] * v1[1]) + P[26] * (v0[2] * v1[5] - v0[5] * v1[2]))) / std::pow(v0[1] * v0[5] - v0[2] * v0[4], 2);
			jacobians[0][29] = (P[24] * (v0[1] * (v0[0] * v0[5] - v0[2] * v0[3] + v0[2] * v1[3] - v0[5] * v1[0]) + (-v0[0] + v1[0]) * (v0[1] * v0[5] - v0[2] * v0[4])) + P[25] * (v0[1] * (v0[2] * v1[4] - v0[5] * v1[1]) + v1[1] * (v0[1] * v0[5] - v0[2] * v0[4])) + P[26] * (v0[1] * (v0[2] * v1[5] - v0[5] * v1[2]) + v1[2] * (v0[1] * v0[5] - v0[2] * v0[4])) - v0[1] * (P[27] * (v0[0] * v0[4] - v0[1] * v0[3] + v0[1] * v1[3] - v0[4] * v1[0]) + P[28] * (v0[1] * v1[4] - v0[4] * v1[1]) + P[29] * (v0[1] * v1[5] - v0[4] * v1[2]))) / std::pow(v0[1] * v0[5] - v0[2] * v0[4], 2);
			jacobians[0][30] = (-P[30] * v0[5] + P[33] * v0[4]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[0][31] = (-P[33] * (v0[5] * (v0[0] * v0[4] - v0[1] * v0[3] + v0[1] * v1[3] - v0[4] * v1[0]) + (v0[3] - v1[3]) * (v0[1] * v0[5] - v0[2] * v0[4])) + P[34] * (v0[5] * (-v0[1] * v1[4] + v0[4] * v1[1]) + v1[4] * (v0[1] * v0[5] - v0[2] * v0[4])) + P[35] * (v0[5] * (-v0[1] * v1[5] + v0[4] * v1[2]) + v1[5] * (v0[1] * v0[5] - v0[2] * v0[4])) + v0[5] * (P[30] * (v0[0] * v0[5] - v0[2] * v0[3] + v0[2] * v1[3] - v0[5] * v1[0]) + P[31] * (v0[2] * v1[4] - v0[5] * v1[1]) + P[32] * (v0[2] * v1[5] - v0[5] * v1[2]))) / std::pow(v0[1] * v0[5] - v0[2] * v0[4], 2);
			jacobians[0][32] = (-P[30] * (v0[4] * (v0[0] * v0[5] - v0[2] * v0[3] + v0[2] * v1[3] - v0[5] * v1[0]) + (-v0[3] + v1[3]) * (v0[1] * v0[5] - v0[2] * v0[4])) - P[31] * (v0[4] * (v0[2] * v1[4] - v0[5] * v1[1]) + v1[4] * (v0[1] * v0[5] - v0[2] * v0[4])) - P[32] * (v0[4] * (v0[2] * v1[5] - v0[5] * v1[2]) + v1[5] * (v0[1] * v0[5] - v0[2] * v0[4])) + v0[4] * (P[33] * (v0[0] * v0[4] - v0[1] * v0[3] + v0[1] * v1[3] - v0[4] * v1[0]) + P[34] * (v0[1] * v1[4] - v0[4] * v1[1]) + P[35] * (v0[1] * v1[5] - v0[4] * v1[2]))) / std::pow(v0[1] * v0[5] - v0[2] * v0[4], 2);
			jacobians[0][33] = (P[30] * v0[2] - P[33] * v0[1]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[0][34] = (P[33] * (v0[2] * (v0[0] * v0[4] - v0[1] * v0[3] + v0[1] * v1[3] - v0[4] * v1[0]) + (v0[0] - v1[0]) * (v0[1] * v0[5] - v0[2] * v0[4])) - P[34] * (v0[2] * (-v0[1] * v1[4] + v0[4] * v1[1]) + v1[1] * (v0[1] * v0[5] - v0[2] * v0[4])) - P[35] * (v0[2] * (-v0[1] * v1[5] + v0[4] * v1[2]) + v1[2] * (v0[1] * v0[5] - v0[2] * v0[4])) - v0[2] * (P[30] * (v0[0] * v0[5] - v0[2] * v0[3] + v0[2] * v1[3] - v0[5] * v1[0]) + P[31] * (v0[2] * v1[4] - v0[5] * v1[1]) + P[32] * (v0[2] * v1[5] - v0[5] * v1[2]))) / std::pow(v0[1] * v0[5] - v0[2] * v0[4], 2);
			jacobians[0][35] = (P[30] * (v0[1] * (v0[0] * v0[5] - v0[2] * v0[3] + v0[2] * v1[3] - v0[5] * v1[0]) + (-v0[0] + v1[0]) * (v0[1] * v0[5] - v0[2] * v0[4])) + P[31] * (v0[1] * (v0[2] * v1[4] - v0[5] * v1[1]) + v1[1] * (v0[1] * v0[5] - v0[2] * v0[4])) + P[32] * (v0[1] * (v0[2] * v1[5] - v0[5] * v1[2]) + v1[2] * (v0[1] * v0[5] - v0[2] * v0[4])) - v0[1] * (P[33] * (v0[0] * v0[4] - v0[1] * v0[3] + v0[1] * v1[3] - v0[4] * v1[0]) + P[34] * (v0[1] * v1[4] - v0[4] * v1[1]) + P[35] * (v0[1] * v1[5] - v0[4] * v1[2]))) / std::pow(v0[1] * v0[5] - v0[2] * v0[4], 2);

		}
		if (jacobians[1]) {
			jacobians[1][0] = (P[0] * v0[5] - P[3] * v0[4]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[1][1] = (P[1] * v0[5] - P[4] * v0[4]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[1][2] = (P[2] * v0[5] - P[5] * v0[4]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[1][3] = (-P[0] * v0[2] + P[3] * v0[1]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[1][4] = (-P[1] * v0[2] + P[4] * v0[1]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[1][5] = (-P[2] * v0[2] + P[5] * v0[1]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[1][6] = (P[6] * v0[5] - P[9] * v0[4]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[1][7] = (-P[10] * v0[4] + P[7] * v0[5]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[1][8] = (-P[11] * v0[4] + P[8] * v0[5]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[1][9] = (-P[6] * v0[2] + P[9] * v0[1]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[1][10] = (P[10] * v0[1] - P[7] * v0[2]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[1][11] = (P[11] * v0[1] - P[8] * v0[2]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[1][12] = (P[12] * v0[5] - P[15] * v0[4]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[1][13] = (P[13] * v0[5] - P[16] * v0[4]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[1][14] = (P[14] * v0[5] - P[17] * v0[4]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[1][15] = (-P[12] * v0[2] + P[15] * v0[1]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[1][16] = (-P[13] * v0[2] + P[16] * v0[1]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[1][17] = (-P[14] * v0[2] + P[17] * v0[1]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[1][18] = (P[18] * v0[5] - P[21] * v0[4]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[1][19] = (P[19] * v0[5] - P[22] * v0[4]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[1][20] = (P[20] * v0[5] - P[23] * v0[4]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[1][21] = (-P[18] * v0[2] + P[21] * v0[1]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[1][22] = (-P[19] * v0[2] + P[22] * v0[1]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[1][23] = (-P[20] * v0[2] + P[23] * v0[1]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[1][24] = (P[24] * v0[5] - P[27] * v0[4]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[1][25] = (P[25] * v0[5] - P[28] * v0[4]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[1][26] = (P[26] * v0[5] - P[29] * v0[4]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[1][27] = (-P[24] * v0[2] + P[27] * v0[1]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[1][28] = (-P[25] * v0[2] + P[28] * v0[1]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[1][29] = (-P[26] * v0[2] + P[29] * v0[1]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[1][30] = (P[30] * v0[5] - P[33] * v0[4]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[1][31] = (P[31] * v0[5] - P[34] * v0[4]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[1][32] = (P[32] * v0[5] - P[35] * v0[4]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[1][33] = (-P[30] * v0[2] + P[33] * v0[1]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[1][34] = (-P[31] * v0[2] + P[34] * v0[1]) / (v0[1] * v0[5] - v0[2] * v0[4]);
			jacobians[1][35] = (-P[32] * v0[2] + P[35] * v0[1]) / (v0[1] * v0[5] - v0[2] * v0[4]);
		}
		return true;
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

			//ceres::CostFunction* cost_function = AdditiveWithPrior_PoseError::Create(e.RelParaVec.data(), e.InvCovMatrix.data());
			ceres::CostFunction* cost_function = new Affine2WithPrior_Additive_Analytic(e.RelParaVec.data(), e.InvCovMatrix.data());
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

/****************************************************************************************************
x2 = Node_para[0] + Node_para[1] * x1 + Node_para[2] * y1;
y2 = Node_para[3] + Node_para[4] * x1 + Node_para[5] * y1;
******************************************************************************************************/

bool MoAve_Affine::_LocalSolver(int index_i)
{

	if (_num_para_per_node != 6)
	{
		std::cout << "number of parameter per node does not match the problem" << std::endl;
		return false;
	}

	std::vector<double> temvec(_num_para_per_node, 0);
	std::vector<double> temWeightVec(_num_para_per_node, 0);

	for (int p1 = 0; p1 < (*_Con_Graph)[index_i].size(); p1++)
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
		XfmComposite(otherNodePara.data(), invRelPara, estRef);

		// 0
		CurTemWeight = (*_Con_Graph)[index_i][p1].InvCovMatrix[0 * _num_para_per_node + 0] * (*_Con_Graph)[index_i][p1].LinkWeight * Multiplier;

		temvec[0] += CurTemWeight * estRef[0];

		temWeightVec[0] += CurTemWeight;

		// 1
		CurTemWeight = (*_Con_Graph)[index_i][p1].InvCovMatrix[1 * _num_para_per_node + 1] * (*_Con_Graph)[index_i][p1].LinkWeight * Multiplier;

		temvec[1] += CurTemWeight * estRef[1];

		temWeightVec[1] += CurTemWeight;

		// 2
		CurTemWeight = (*_Con_Graph)[index_i][p1].InvCovMatrix[2 * _num_para_per_node + 2] * (*_Con_Graph)[index_i][p1].LinkWeight * Multiplier;

		temvec[2] += CurTemWeight * estRef[2];

		temWeightVec[2] += CurTemWeight;


		// 3
		CurTemWeight = (*_Con_Graph)[index_i][p1].InvCovMatrix[3 * _num_para_per_node + 3] * (*_Con_Graph)[index_i][p1].LinkWeight * Multiplier;

		temvec[3] += CurTemWeight * estRef[3];

		temWeightVec[3] += CurTemWeight;

		// 4
		CurTemWeight = (*_Con_Graph)[index_i][p1].InvCovMatrix[4 * _num_para_per_node + 4] * (*_Con_Graph)[index_i][p1].LinkWeight * Multiplier;

		temvec[4] += CurTemWeight * estRef[4];

		temWeightVec[4] += CurTemWeight;

		// 5
		CurTemWeight = (*_Con_Graph)[index_i][p1].InvCovMatrix[5 * _num_para_per_node + 5] * (*_Con_Graph)[index_i][p1].LinkWeight * Multiplier;

		temvec[5] += CurTemWeight * estRef[5];

		temWeightVec[5] += CurTemWeight;

	}// end for p1

	for (int p2 = 0; p2 < _num_para_per_node; p2++)
		_Node_Para_list_buf[index_i][p2] = temvec[p2] / temWeightVec[p2];             // averaging the current updates

	return true;
}