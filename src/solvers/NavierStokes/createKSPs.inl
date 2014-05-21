template <PetscInt dim>
void NavierStokesSolver<dim>::createKSPs()
{
	PetscErrorCode ierr;
	ierr = KSPCreate(PETSC_COMM_WORLD, &ksp1); CHKERRV(ierr);
	ierr = KSPSetTolerances(ksp1, 1e-6, 0.0, PETSC_DEFAULT, PETSC_DEFAULT); CHKERRV(ierr);
	ierr = KSPSetOperators(ksp1, A, A, DIFFERENT_NONZERO_PATTERN); CHKERRV(ierr);
	ierr = KSPSetFromOptions(ksp1); CHKERRV(ierr);
}