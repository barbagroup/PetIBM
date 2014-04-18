#include "NavierStokesSolver.h"
#include <iostream>
#include <string>

template <int dim>
NavierStokesSolver<dim>* NavierStokesSolver<dim>::createSolver(SimulationParameters &SP, CartesianMesh<dim> &CM)
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
	solver->simParams = &SP;
	solver->mesh      = &CM;
	
	std::cout << "\nSolver type selected: " << solver->name() << "\n" << std::endl;
	
	return solver;
}

template class NavierStokesSolver<2>;
template class NavierStokesSolver<3>;
