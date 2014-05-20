#include <createSolver.h>

template <PetscInt dim>
std::unique_ptr< NavierStokesSolver<dim> > createSolver(FlowDescription *FD, SimulationParameters *SP, CartesianMesh *CM)
{
	PetscPrintf(PETSC_COMM_WORLD, "gamma: %f, zeta: %f, alphaExplicit: %f, alphaImplicit: %f\n", SP->gamma, SP->zeta, SP->alphaExplicit, SP->alphaImplicit);
	PetscPrintf(PETSC_COMM_WORLD, "Solver type selected: ");
	switch(SP->solverType)
	{
		case NAVIER_STOKES:
			PetscPrintf(PETSC_COMM_WORLD, "Navier-Stokes\n");
			return std::unique_ptr< NavierStokesSolver<dim> >(new NavierStokesSolver<dim>(FD, SP, CM));
			break;
		default:
			PetscPrintf(PETSC_COMM_WORLD, "Unrecognised solver!\n");
			return NULL;
	}
}

template std::unique_ptr< NavierStokesSolver<2> > createSolver(FlowDescription *FD, SimulationParameters *SP, CartesianMesh *CM);
template std::unique_ptr< NavierStokesSolver<3> > createSolver(FlowDescription *FD, SimulationParameters *SP, CartesianMesh *CM);