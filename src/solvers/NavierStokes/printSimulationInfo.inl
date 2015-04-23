/***************************************************************************//**
 * \file printSimulationInfo.inl
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the method `printSimulationInfo` 
 *        of the class \c NavierStokesSolver.
 */

/**
 * \brief Prints the simulation parameters.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::printSimulationInfo()
{
	PetscErrorCode ierr;
	PetscInt       rank, maxits;
	PC             pc;
	PCType         pcType;
	KSPType        kspType;
	PetscReal      rtol, abstol;
	
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n---------------------------------------\n"); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "Flow\n"); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "---------------------------------------\n"); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "solver: "); CHKERRQ(ierr);
	switch (simParams->solverType)
	{
		case NAVIER_STOKES :
			ierr = PetscPrintf(PETSC_COMM_WORLD, "Navier-Stokes\n"); CHKERRQ(ierr); 
			break;
		case TAIRA_COLONIUS: 
			ierr = PetscPrintf(PETSC_COMM_WORLD, "Taira & Colonius (2007)\n"); CHKERRQ(ierr);
			break;
		default:
			ierr = PetscPrintf(PETSC_COMM_WORLD, "Unrecognized solver!\n"); CHKERRQ(ierr);
			break;
	}
	ierr = PetscPrintf(PETSC_COMM_WORLD, "viscosity: %g\n", flowDesc->nu); CHKERRQ(ierr);
	
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n---------------------------------------\n"); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "Mesh\n"); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "---------------------------------------\n"); CHKERRQ(ierr);
	if (dim == 3)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "size: %d x %d x %d\n", mesh->nx, mesh->ny, mesh->nz); CHKERRQ(ierr);
	}
	else
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "size: %d x %d\n", mesh->nx, mesh->ny); CHKERRQ(ierr);
	}
	
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n---------------------------------------\n"); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "Time-stepping\n"); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "---------------------------------------\n"); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "convection: "); CHKERRQ(ierr);
	switch (simParams->convectionScheme)
	{
		case EULER_EXPLICIT   :
			ierr = PetscPrintf(PETSC_COMM_WORLD, "Explicit Euler\n"); CHKERRQ(ierr);
			break;
		case ADAMS_BASHFORTH_2:
			ierr = PetscPrintf(PETSC_COMM_WORLD, "2nd-order Adams-Bashforth\n"); CHKERRQ(ierr);
			break;
		default:
			break;
	}
	ierr = PetscPrintf(PETSC_COMM_WORLD, "diffusion : "); CHKERRQ(ierr);
	switch (simParams->diffusionScheme)
	{
		case EULER_EXPLICIT:
			ierr = PetscPrintf(PETSC_COMM_WORLD, "Explicit Euler\n"); CHKERRQ(ierr);
			break;
		case EULER_IMPLICIT:
			ierr = PetscPrintf(PETSC_COMM_WORLD, "Implicit Euler\n"); CHKERRQ(ierr);
			break;
		case CRANK_NICOLSON:
			ierr = PetscPrintf(PETSC_COMM_WORLD, "Crank-Nicolson\n"); CHKERRQ(ierr);
			break;
		default:
			break;
	}
	ierr = PetscPrintf(PETSC_COMM_WORLD, "time-increment      : %g\n", simParams->dt); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "starting time-step  : %d\n", simParams->startStep); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "number of time-steps: %d\n", simParams->nt); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "saving-interval     : %d\n", simParams->nsave); CHKERRQ(ierr);

	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n---------------------------------------\n"); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "Linear system for intermediate velocity\n"); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "---------------------------------------\n"); CHKERRQ(ierr);
	ierr = KSPGetType(ksp1, &kspType); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "solver: %s\n", kspType); CHKERRQ(ierr);
	ierr = KSPGetPC(ksp1, &pc); CHKERRQ(ierr);
	ierr = PCGetType(pc, &pcType); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "preconditioner: %s\n", pcType); CHKERRQ(ierr);
	ierr = KSPGetTolerances(ksp1, &rtol, &abstol, NULL, &maxits); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "relative tolerance: %g\n", rtol); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "absolute tolerance: %g\n", abstol); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "maximum iterations: %d\n", maxits); CHKERRQ(ierr);

	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n---------------------------------------\n"); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "Linear system for pressure-force\n"); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "---------------------------------------\n"); CHKERRQ(ierr);
	ierr = KSPGetType(ksp2, &kspType); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "solver: %s\n", kspType); CHKERRQ(ierr);
	ierr = KSPGetPC(ksp2, &pc); CHKERRQ(ierr);
	ierr = PCGetType(pc, &pcType); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "preconditioner: %s\n", pcType); CHKERRQ(ierr);
	ierr = KSPGetTolerances(ksp2, &rtol, &abstol, NULL, &maxits); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "relative tolerance: %g\n", rtol); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "absolute tolerance: %g\n", abstol); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "maximum iterations: %d\n", maxits); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n---------------------------------------\n"); CHKERRQ(ierr);

	return 0;
}