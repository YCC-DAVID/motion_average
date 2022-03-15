#pragma once;
#include "Motion_Average.h"




/*******************************************************************************************
two dimentional affine transformation
A1X1 - A2X2 = 0; ->
inv(A2)A1X1 - X2 = 0;  A21 = inv(A2)A1;
A1 = A2 * A21;

ConG_Unit.RelParaVec stores A21
A1 refers to reference
A2 refers to otherID

********************************************************************************************/


class MoAve_Affine : public MotionAverage_T
{

public:

	MoAve_Affine();
	~MoAve_Affine(); 

	bool DirectSolver_Run();

private:
	bool _LocalSolver(int index_i);

};