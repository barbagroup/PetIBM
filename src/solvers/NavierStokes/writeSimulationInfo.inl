template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::writeSimulationInfo()
{
	PetscErrorCode ierr;
	PetscInt       rank, maxits;
	PC             pc;
	PCType         pcType;
	KSPType        kspType;
	PetscReal      rtol;
	
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n---------------------------------------\n"); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "Flow\n"); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "---------------------------------------\n"); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "Solver type: "); CHKERRQ(ierr);
	switch(simParams->solverType)
	{
		case NAVIER_STOKES : ierr = PetscPrintf(PETSC_COMM_WORLD, "Navier-Stokes\n"); CHKERRQ(ierr); break;
		case TAIRA_COLONIUS: ierr = PetscPrintf(PETSC_COMM_WORLD, "Taira & Colonius (2007)\n"); CHKERRQ(ierr); break;
		default: ierr = PetscPrintf(PETSC_COMM_WORLD, "Unrecognized solver!\n"); CHKERRQ(ierr); break;
	}
	ierr = PetscPrintf(PETSC_COMM_WORLD, "nu: %f\n", flowDesc->nu); CHKERRQ(ierr);
	
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n---------------------------------------\n"); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "Mesh\n"); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "---------------------------------------\n"); CHKERRQ(ierr);
	if(dim == 3)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "Size: %d x %d x %d\n", mesh->nx, mesh->ny, mesh->nz); CHKERRQ(ierr);
	}
	else
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "Size: %d x %d\n", mesh->nx, mesh->ny); CHKERRQ(ierr);
	}
	
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n---------------------------------------\n"); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "Time-stepping schemes\n"); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "---------------------------------------\n"); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "Convection: "); CHKERRQ(ierr);
	switch(simParams->convectionScheme)
	{
		case EULER_EXPLICIT   : ierr = PetscPrintf(PETSC_COMM_WORLD, "Explicit Euler\n"); CHKERRQ(ierr); break;
		case ADAMS_BASHFORTH_2: ierr = PetscPrintf(PETSC_COMM_WORLD, "2nd-order Adams-Bashforth\n"); CHKERRQ(ierr); break;
		default: break;
	}
	ierr = PetscPrintf(PETSC_COMM_WORLD, "Diffusion : "); CHKERRQ(ierr);
	switch(simParams->diffusionScheme)
	{
		case EULER_EXPLICIT: ierr = PetscPrintf(PETSC_COMM_WORLD, "Explicit Euler\n"); CHKERRQ(ierr); break;
		case EULER_IMPLICIT: ierr = PetscPrintf(PETSC_COMM_WORLD, "Implicit Euler\n"); CHKERRQ(ierr); break;
		case CRANK_NICOLSON: ierr = PetscPrintf(PETSC_COMM_WORLD, "Crank-Nicolson\n"); CHKERRQ(ierr); break;
		default: break;
	}
	ierr = PetscPrintf(PETSC_COMM_WORLD, "startStep: %d\n", simParams->startStep); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "nt       : %d\n", simParams->nt); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "nsave    : %d\n", simParams->nsave); CHKERRQ(ierr);

	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n---------------------------------------\n"); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "Linear system for intermediate velocity\n"); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "---------------------------------------\n"); CHKERRQ(ierr);
	ierr = KSPGetType(ksp1, &kspType); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "Solver: %s\n", kspType); CHKERRQ(ierr);
	ierr = KSPGetPC(ksp1, &pc); CHKERRQ(ierr);
	ierr = PCGetType(pc, &pcType); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "Preconditioner: %s\n", pcType); CHKERRQ(ierr);
	ierr = KSPGetTolerances(ksp1, &rtol, NULL, NULL, &maxits); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "Relative tolerance: %f\n", rtol); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "Maximum iterations: %d\n", maxits); CHKERRQ(ierr);

	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n---------------------------------------\n"); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "Linear system for Poisson step\n"); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "---------------------------------------\n"); CHKERRQ(ierr);
	ierr = KSPGetType(ksp2, &kspType); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "Solver: %s\n", kspType); CHKERRQ(ierr);
	ierr = KSPGetPC(ksp2, &pc); CHKERRQ(ierr);
	ierr = PCGetType(pc, &pcType); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "Preconditioner: %s\n", pcType); CHKERRQ(ierr);
	ierr = KSPGetTolerances(ksp2, &rtol, NULL, NULL, &maxits); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "Relative tolerance: %f\n", rtol); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "Maximum iterations: %d\n", maxits); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "\n---------------------------------------\n"); CHKERRQ(ierr);

	if(rank==0)
	{
		std::ofstream f(caseFolder+"/simulationInfo.txt");
		f << "-nx\t" << mesh->nx << '\n';
		f << "-ny\t" << mesh->ny << '\n';
		if(dim == 3)
		{
			f << "-nz\t" << mesh->nz << '\n';
		}
		f << "-startStep\t" << simParams->startStep << '\n';
		f << "-nt\t" << simParams->nt << '\n';
		f << "-nsave\t" << simParams->nsave << '\n';
		f << "-dt\t" << simParams->dt << '\n';
		(flowDesc->bc[0][XPLUS].type==PERIODIC)? f << "-xperiodic\tTrue\n" : f << "-xperiodic\tFalse\n";
		(flowDesc->bc[0][YPLUS].type==PERIODIC)? f << "-yperiodic\tTrue\n" : f << "-yperiodic\tFalse\n";
		if(dim == 3)
		{
			(flowDesc->bc[0][ZPLUS].type==PERIODIC)? f << "-zperiodic\tTrue\n" : f << "-zperiodic\tFalse\n";
		}
		f.close();
	}

	return 0;
}