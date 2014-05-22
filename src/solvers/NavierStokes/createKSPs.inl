template <PetscInt dim>
void NavierStokesSolver<dim>::createKSPs()
{
	PetscErrorCode ierr;
	
	// linear system for the intermediate velocity
	ierr = KSPCreate(PETSC_COMM_WORLD, &ksp1); CHKERRV(ierr);
	ierr = KSPSetTolerances(ksp1, 1e-6, 0.0, PETSC_DEFAULT, PETSC_DEFAULT); CHKERRV(ierr);
	ierr = KSPSetOperators(ksp1, A, A, DIFFERENT_NONZERO_PATTERN); CHKERRV(ierr);
	ierr = KSPSetInitialGuessNonzero(ksp1, PETSC_TRUE); CHKERRV(ierr);
	ierr = KSPSetFromOptions(ksp1); CHKERRV(ierr);

	// linear system for the Poisson solver
	ierr = KSPCreate(PETSC_COMM_WORLD, &ksp2); CHKERRV(ierr);
	ierr = KSPSetTolerances(ksp2, 1e-6, 0.0, PETSC_DEFAULT, PETSC_DEFAULT); CHKERRV(ierr);
	ierr = KSPSetOperators(ksp2, QTBNQ, QTBNQ, DIFFERENT_NONZERO_PATTERN); CHKERRV(ierr);
	ierr = KSPSetInitialGuessNonzero(ksp2, PETSC_TRUE); CHKERRV(ierr);
	ierr = KSPSetFromOptions(ksp2); CHKERRV(ierr);
}