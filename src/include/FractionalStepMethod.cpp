#include "FractionalStepMethod.h"

void FractionalStepMethod::solveIntermediateVelocity()
{
	PetscErrorCode ierr;
	ierr = KSPSolve(ksp1, rhs1, qStar); CHKERRV(ierr);
}

void FractionalStepMethod::solvePoissonSystem()
{
	PetscErrorCode ierr;
	ierr = KSPSolve(ksp2, rhs2, phi); CHKERRV(ierr);
}

void FractionalStepMethod::projectionStep()
{
	
}
