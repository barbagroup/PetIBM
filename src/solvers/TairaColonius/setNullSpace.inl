template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::setNullSpace()
{
	PetscErrorCode ierr;
	MatNullSpace   nsp;
	Vec            phiPortion;
	PetscInt       numPhi;

	ierr = VecSet(nullSpaceVec, 0.0); CHKERRQ(ierr);
	ierr = DMCompositeGetAccess(NavierStokesSolver<dim>::lambdaPack, nullSpaceVec,  &phiPortion, NULL); CHKERRQ(ierr);
	ierr = VecGetSize(phiPortion, &numPhi);
	ierr = VecSet(phiPortion, 1.0/sqrt(numPhi)); CHKERRQ(ierr);
	ierr = DMCompositeRestoreAccess(NavierStokesSolver<dim>::lambdaPack, nullSpaceVec,  &phiPortion, NULL); CHKERRQ(ierr);
	
	ierr = MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_FALSE, 1, &nullSpaceVec, &nsp); CHKERRQ(ierr);
	ierr = KSPSetNullSpace(NavierStokesSolver<dim>::ksp2, nsp);
	ierr = MatNullSpaceDestroy(&nsp);

	return 0;
}