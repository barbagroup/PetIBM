#include "TairaColoniusSolver.h"
#include <petscsys.h>
#include <petscdmcomposite.h>
#include <iostream>
#include <string>

#include "TairaColonius/createDMs.inl"
#include "TairaColonius/generateBNQ.inl"
#include "TairaColonius/generateET.inl"
#include "TairaColonius/generateR2.inl"
#include "TairaColonius/initialiseBodies.inl"
#include "TairaColonius/createGlobalMappingBodies.inl"

template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::initialise()
{
	PetscErrorCode ierr;

	initialiseBodies();
	ierr = createDMs(); CHKERRQ(ierr);
	ierr = NavierStokesSolver<dim>::createVecs(); CHKERRQ(ierr);

	NavierStokesSolver<dim>::initialiseMeshSpacings();
	ierr = NavierStokesSolver<dim>::initialiseFluxes(); CHKERRQ(ierr);
	ierr = NavierStokesSolver<dim>::updateBoundaryGhosts(); CHKERRQ(ierr);

	ierr = createGlobalMappingBodies(); CHKERRQ(ierr);
	ierr = NavierStokesSolver<dim>::createLocalToGlobalMappingsFluxes(); CHKERRQ(ierr);
	ierr = NavierStokesSolver<dim>::createLocalToGlobalMappingsLambda(); CHKERRQ(ierr);

	ierr = NavierStokesSolver<dim>::generateDiagonalMatrices(); CHKERRQ(ierr);
	ierr = NavierStokesSolver<dim>::generateA(); CHKERRQ(ierr);
	ierr = generateBNQ(); CHKERRQ(ierr);
	ierr = NavierStokesSolver<dim>::generateQTBNQ(); CHKERRQ(ierr);
	ierr = NavierStokesSolver<dim>::createKSPs(); CHKERRQ(ierr);
	ierr = generateET(); CHKERRQ(ierr);

	return 0;
}

template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::finalise()
{
	PetscErrorCode ierr;

	ierr = NavierStokesSolver<dim>::finalise();

	// DMs
	if(bda!=PETSC_NULL) {ierr = DMDestroy(&bda); CHKERRQ(ierr);}
	// Mats
	if(ET!=PETSC_NULL)  {ierr = MatDestroy(&ET); CHKERRQ(ierr);}

	return 0;
}

template <PetscInt dim>
PetscReal TairaColoniusSolver<dim>::dhRoma(PetscReal x, PetscReal h)
{
	PetscReal r = fabs(x)/h;
	if(r>1.5)
		return 0.0;
	else if(r>0.5 && r<=1.5)
		return 1.0/(6*h)*( 5.0 - 3.0*r - sqrt(-3.0*(1-r)*(1-r) + 1.0) );
	else
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

template class TairaColoniusSolver<2>;
template class TairaColoniusSolver<3>;