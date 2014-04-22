#include "NavierStokesSolver.h"
#include <petscsys.h>
#include <iostream>
#include <string>

// this factory pattern causes a memory leak. Try to fix it later
template <PetscInt dim>
NavierStokesSolver<dim>* NavierStokesSolver<dim>::createSolver(FlowDescription &FD, SimulationParameters &SP, CartesianMesh &CM)
{
	NavierStokesSolver<dim> *solver = NULL;
	switch(SP.solverType)
	{
		case NAVIER_STOKES:
			solver = new NavierStokesSolver<dim>;
			break;
		default:
			std::cout << "Unrecognised solver!\n";
	}
	solver->flowDesc  = &FD;
	solver->simParams = &SP;
	solver->mesh      = &CM;
	
	PetscPrintf(PETSC_COMM_WORLD, "Solver type selected: %s\n", solver->name().c_str());
	
	return solver;
}

template <PetscInt dim>
void NavierStokesSolver<dim>::initialise()
{
	// Initialise fluxes
	fluxVecsCreate();
	fluxVecsInitialise();
	updateBoundaryGhosts();
}

#include "NavierStokes/fluxVecsCreate.inl"
#include "NavierStokes/fluxVecsInitialise.inl"
#include "NavierStokes/updateBoundaryGhosts.inl"

template class NavierStokesSolver<2>;
template class NavierStokesSolver<3>;
