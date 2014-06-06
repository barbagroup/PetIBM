template <PetscInt dim>
void NavierStokesSolver<dim>::createKSPs()
{
	PetscErrorCode ierr;
	
	// linear system for the intermediate velocity
	ierr = KSPCreate(PETSC_COMM_WORLD, &ksp1); CHKERRV(ierr);
	ierr = KSPSetTolerances(ksp1, 1e-6, 0.0, PETSC_DEFAULT, PETSC_DEFAULT); CHKERRV(ierr);
	ierr = KSPSetOperators(ksp1, A, A, SAME_NONZERO_PATTERN); CHKERRV(ierr);
	ierr = KSPSetInitialGuessNonzero(ksp1, PETSC_TRUE); CHKERRV(ierr);
	ierr = KSPSetFromOptions(ksp1); CHKERRV(ierr);

	// linear system for the Poisson solver
	ierr = KSPCreate(PETSC_COMM_WORLD, &ksp2); CHKERRV(ierr);
	ierr = KSPSetTolerances(ksp2, 1e-6, 0.0, PETSC_DEFAULT, PETSC_DEFAULT); CHKERRV(ierr);
	ierr = KSPSetOperators(ksp2, QTBNQ, QTBNQ, SAME_NONZERO_PATTERN); CHKERRV(ierr);
	ierr = KSPSetInitialGuessNonzero(ksp2, PETSC_TRUE); CHKERRV(ierr);
	ierr = KSPGetPC(ksp2, &pc2); CHKERRV(ierr);
	ierr = KSPSetType(ksp2, KSPCG); CHKERRV(ierr);
	ierr = PCSetType(pc2, PCGAMG); CHKERRV(ierr);
	ierr = KSPSetFromOptions(ksp2); CHKERRV(ierr);
}