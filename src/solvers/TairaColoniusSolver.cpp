#include "TairaColoniusSolver.h"
#include <petscsys.h>
#include <petscdmcomposite.h>
#include <iostream>
#include <string>

#include "TairaColonius/createDMs.inl"
#include "TairaColonius/generateBNQ.inl"
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
	ierr = createGlobalMappingBodies(); CHKERRQ(ierr);
	ierr = NavierStokesSolver<dim>::createLocalToGlobalMappingsFluxes(); CHKERRQ(ierr);
	ierr = NavierStokesSolver<dim>::createLocalToGlobalMappingsLambda(); CHKERRQ(ierr);
	NavierStokesSolver<dim>::initialiseMeshSpacings();
	ierr = NavierStokesSolver<dim>::initialiseFluxes(); CHKERRQ(ierr);
	ierr = NavierStokesSolver<dim>::updateBoundaryGhosts(); CHKERRQ(ierr);

	ierr = NavierStokesSolver<dim>::generateDiagonalMatrices(); CHKERRQ(ierr);
	ierr = NavierStokesSolver<dim>::generateA(); CHKERRQ(ierr);
	ierr = generateBNQ(); CHKERRQ(ierr);
	ierr = NavierStokesSolver<dim>::generateQTBNQ(); CHKERRQ(ierr);
	ierr = NavierStokesSolver<dim>::createKSPs(); CHKERRQ(ierr);

	return 0;
}

template class TairaColoniusSolver<2>;
template class TairaColoniusSolver<3>;