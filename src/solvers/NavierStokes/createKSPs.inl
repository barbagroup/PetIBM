template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::createKSPs()
{
	PetscErrorCode ierr;
	
	// linear system for the intermediate velocity
	ierr = KSPCreate(PETSC_COMM_WORLD, &ksp1); CHKERRQ(ierr);
	ierr = KSPSetOptionsPrefix(ksp1, "sys1_"); CHKERRQ(ierr);
	ierr = KSPSetTolerances(ksp1, 1e-6, 0.0, PETSC_DEFAULT, PETSC_DEFAULT); CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp1, A, A, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
	ierr = KSPSetInitialGuessNonzero(ksp1, PETSC_TRUE); CHKERRQ(ierr);
	ierr = KSPSetType(ksp1, KSPCG); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp1); CHKERRQ(ierr);

	// linear system for the Poisson solver
	ierr = KSPCreate(PETSC_COMM_WORLD, &ksp2); CHKERRQ(ierr);
	ierr = KSPSetOptionsPrefix(ksp2, "sys2_"); CHKERRQ(ierr);
	ierr = KSPSetTolerances(ksp2, 1e-6, 0.0, PETSC_DEFAULT, PETSC_DEFAULT); CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp2, QTBNQ, QTBNQ, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
	ierr = KSPSetInitialGuessNonzero(ksp2, PETSC_TRUE); CHKERRQ(ierr);
	ierr = KSPSetType(ksp2, KSPCG); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp2); CHKERRQ(ierr);

	return 0;
}