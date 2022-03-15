#pragma once;
#include "Motion_Average.h"





/************************************************************************************************
This is the implementation of motion average for Zero order problem defined as follows
P_i = {Px_0,Px_1,Px_2,Px_3,...,Px_n}, 
     where F_ij(P_i,P_j,P_ij,cov(x_i,x_j)) = 0 
  => featureless constraints
     P_i - P_j = P_ij, cov defines the contribution of each element of P_i towards the error of this constraint
	 since this is linear, P_i = P_ij + P_j gives zero error, this cov will not be used here, 
	 but in the aggregated InvCov to get Mahalanobis distance
*************************************************************************************************/
class MoAve_0 : public MotionAverage_T
{

public:
	MoAve_0();
	~MoAve_0();


	bool DirectSolver_Run();

private: 
	bool _LocalSolver(int index_i);
	
};