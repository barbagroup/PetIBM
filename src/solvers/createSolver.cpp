#include "createSolver.h"
#include "TairaColoniusSolver.h"

template <PetscInt dim>
std::unique_ptr< NavierStokesSolver<dim> > createSolver(std::string folder, FlowDescription *FD, SimulationParameters *SP, CartesianMesh *CM)
{
	if(FD->dimensions!=dim)
	{
		PetscPrintf(PETSC_COMM_WORLD, "ERROR: Check the number of dimensions in the input files.\n");
		exit(0);
	}
	PetscPrintf(PETSC_COMM_WORLD, "gamma: %f, zeta: %f, alphaExplicit: %f, alphaImplicit: %f\n", SP->gamma, SP->zeta, SP->alphaExplicit, SP->alphaImplicit);
	PetscPrintf(PETSC_COMM_WORLD, "Solver type selected: ");
	switch(SP->solverType)
	{
		case NAVIER_STOKES:
			PetscPrintf(PETSC_COMM_WORLD, "Navier-Stokes\n");
			return std::unique_ptr< NavierStokesSolver<dim> >(new NavierStokesSolver<dim>(folder, FD, SP, CM));
			break;
		case TAIRA_COLONIUS:
			PetscPrintf(PETSC_COMM_WORLD, "Taira and Colonius\n");
			return std::unique_ptr< TairaColoniusSolver<dim> >(new TairaColoniusSolver<dim>(folder, FD, SP, CM));
			break;
		default:
			PetscPrintf(PETSC_COMM_WORLD, "Unrecognised solver!\n");
			return NULL;
	}
}

template std::unique_ptr< NavierStokesSolver<2> > createSolver(std::string folder, FlowDescription *FD, SimulationParameters *SP, CartesianMesh *CM);
template std::unique_ptr< NavierStokesSolver<3> > createSolver(std::string folder, FlowDescription *FD, SimulationParameters *SP, CartesianMesh *CM);