#include "TairaColoniusSolver.h"
#include "yaml-cpp/yaml.h"
#include <petscdmcomposite.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <sys/stat.h>

template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::initialize()
{
	PetscErrorCode ierr;

	ierr = PetscLogStagePush(NavierStokesSolver<dim>::stageInitialize); CHKERRQ(ierr);
	ierr = initializeBodies(); CHKERRQ(ierr);
	ierr = calculateCellIndices(); CHKERRQ(ierr);
	ierr = createDMs(); CHKERRQ(ierr);
	ierr = createGlobalMappingBodies(); CHKERRQ(ierr);
	ierr = NavierStokesSolver<dim>::initializeCommon(); CHKERRQ(ierr);
	ierr = VecDuplicate(NavierStokesSolver<dim>::q, &temp); CHKERRQ(ierr);
	ierr = PetscLogStagePop(); CHKERRQ(ierr);

	return 0;
}

template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::finalize()
{
	PetscErrorCode ierr;

	ierr = NavierStokesSolver<dim>::finalize();

	// DMs
	if(bda!=PETSC_NULL) {ierr = DMDestroy(&bda); CHKERRQ(ierr);}
	// Mats
	if(ET!=PETSC_NULL)  {ierr = MatDestroy(&ET); CHKERRQ(ierr);}
	// Vecs
	if(temp!=PETSC_NULL){ierr = VecDestroy(&temp); CHKERRQ(ierr);}

	return 0;
}

template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::createDMs()
{
	PetscErrorCode ierr;
	ierr = NavierStokesSolver<dim>::createDMs(); CHKERRQ(ierr);	
	ierr = generateBodyInfo(); CHKERRQ(ierr);
	ierr = DMDACreate1d(PETSC_COMM_WORLD, DMDA_BOUNDARY_NONE, x.size(), dim, 0, &numBoundaryPointsOnProcess.front(), &bda); CHKERRQ(ierr);
	ierr = DMCompositeAddDM(NavierStokesSolver<dim>::lambdaPack, bda); CHKERRQ(ierr);

	return 0;
}

template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::writeData()
{
	NavierStokesSolver<dim>::writeData();
	calculateForce();
	writeForces();

	return 0;
}

template <PetscInt dim>
PetscReal TairaColoniusSolver<dim>::dhRoma(PetscReal x, PetscReal h)
{
	PetscReal r = fabs(x)/h;
	if(r>1.5) return 0.0;
	if(r>0.5 && r<=1.5) return 1.0/(6*h)*( 5.0 - 3.0*r - sqrt(-3.0*(1-r)*(1-r) + 1.0) );
	return 1.0/(3*h)*( 1.0 + sqrt(-3.0*r*r + 1.0) );
}

template <PetscInt dim>
PetscReal TairaColoniusSolver<dim>::delta(PetscReal x, PetscReal y, PetscReal h)
{
	return dhRoma(x, h) * dhRoma(y, h);
}

template <PetscInt dim>
PetscReal TairaColoniusSolver<dim>::delta(PetscReal x, PetscReal y, PetscReal z, PetscReal h)
{
	return dhRoma(x, h) * dhRoma(y, h) * dhRoma(z, h);
}

#include "TairaColonius/calculateCellIndices.inl"
#include "TairaColonius/initializeLambda.inl"
#include "TairaColonius/generateBodyInfo.inl"
#include "TairaColonius/generateBNQ.inl"
#include "TairaColonius/generateR2.inl"
#include "TairaColonius/initializeBodies.inl"
#include "TairaColonius/createGlobalMappingBodies.inl"
#include "TairaColonius/isInfluenced.inl"
#include "TairaColonius/writeLambda.inl"
#include "TairaColonius/calculateForce.inl"
#include "TairaColonius/writeForces.inl"

template class TairaColoniusSolver<2>;
template class TairaColoniusSolver<3>;
