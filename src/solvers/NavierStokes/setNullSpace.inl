template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::setNullSpace()
{
	PetscErrorCode ierr;
	MatNullSpace   nsp;

	ierr = MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, NULL, &nsp); CHKERRQ(ierr);
	ierr = KSPSetNullSpace(ksp2, nsp);
	ierr = MatNullSpaceDestroy(&nsp);

	return 0;
}