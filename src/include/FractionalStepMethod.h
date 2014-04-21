#if !defined(FRACTIONAL_STEP_METHOD_H)
#define FRACTIONAL_STEP_METHOD_H

#include <petscksp.h>

typedef struct
{
	Mat A;
	Mat QT, BNQ;
	Mat QTBNQ;
	Vec BN;
	Vec rhs1, rhs2;
	Vec q, qStar, phi;
	KSP ksp1, ksp2;
	
	void solveIntermediateVelocity();
	void solvePoissonSystem();
	void projectionStep();
public:
	void stepTime();
}FractionalStepMethod;

#endif
