#include "NavierStokesSolver.h"
#include <petscsys.h>
#include <iostream>
#include <string>

// this factory pattern causes a memory leak. Try to fix it later
template <PetscInt dim>
NavierStokesSolver<dim>* NavierStokesSolver<dim>::createSolver(FlowDescription *FD, SimulationParameters *SP, CartesianMesh *CM)
{
	NavierStokesSolver<dim> *solver = NULL;
	switch(SP->solverType)
	{
		case NAVIER_STOKES:
			solver = new NavierStokesSolver<dim>(FD, SP, CM);
			break;
		default:
			std::cout << "Unrecognised solver!\n";
	}
	
	PetscPrintf(PETSC_COMM_WORLD, "Solver type selected: %s\n", solver->name().c_str());
	PetscPrintf(PETSC_COMM_WORLD, "gamma: %f, zeta: %f, alphaExplicit: %f, alphaImplicit: %f\n", solver->simParams->gamma, solver->simParams->zeta, solver->simParams->alphaExplicit, solver->simParams->alphaImplicit);
	
	return solver;
}

template <PetscInt dim>
void NavierStokesSolver<dim>::initialise()
{
	// Initialise fluxes
	createDMs();
	createVecs();
	initialiseFluxes();
	updateBoundaryGhosts();
}

template <PetscInt dim>
void NavierStokesSolver<dim>::stepTime()
{
	timeStep++;
}

template <PetscInt dim>
void NavierStokesSolver<dim>::writeData()
{
}

template <PetscInt dim>
bool NavierStokesSolver<dim>::finished()
{
	return (timeStep < simParams->nt)? false : true;
}

#include "NavierStokes/createDMs.inl"
#include "NavierStokes/createVecs.inl"
#include "NavierStokes/initialiseFluxes.inl"
#include "NavierStokes/updateBoundaryGhosts.inl"
#include "NavierStokes/generateMAndRInv.inl"
#include "NavierStokes/calculateExplicitTerms.inl"

template class NavierStokesSolver<2>;
template class NavierStokesSolver<3>;
