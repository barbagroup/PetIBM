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
void TairaColoniusSolver<dim>::initialise()
{
	initialiseBodies();
	createDMs();
	NavierStokesSolver<dim>::createVecs();
	createGlobalMappingBodies();
	NavierStokesSolver<dim>::createLocalToGlobalMappingsFluxes();
	NavierStokesSolver<dim>::createLocalToGlobalMappingsLambda();
	NavierStokesSolver<dim>::initialiseMeshSpacings();
	NavierStokesSolver<dim>::initialiseFluxes();
	NavierStokesSolver<dim>::updateBoundaryGhosts();

	NavierStokesSolver<dim>::generateDiagonalMatrices();
	NavierStokesSolver<dim>::generateA();
	generateBNQ();
	NavierStokesSolver<dim>::generateQTBNQ();
	NavierStokesSolver<dim>::createKSPs();
}

template class TairaColoniusSolver<2>;
template class TairaColoniusSolver<3>;