/***************************************************************************//**
 * \file setNullSpace.inl
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the method `setNullSpace` of the class `TairaColoniusSolver`.
 */


/**
 * \brief Sets the nullspace of the modified Poisson system.
 */
template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::setNullSpace()
{
  PetscErrorCode ierr;

  ierr = VecSet(nullSpaceVec, 0.0); CHKERRQ(ierr);
  Vec phiPortion;
  ierr = DMCompositeGetAccess(NavierStokesSolver<dim>::lambdaPack, nullSpaceVec,  &phiPortion, NULL); CHKERRQ(ierr);
  PetscInt numPhi;
  ierr = VecGetSize(phiPortion, &numPhi);
  ierr = VecSet(phiPortion, 1.0/sqrt(numPhi)); CHKERRQ(ierr);
  ierr = DMCompositeRestoreAccess(NavierStokesSolver<dim>::lambdaPack, nullSpaceVec,  &phiPortion, NULL); CHKERRQ(ierr);
  MatNullSpace nsp;
  ierr = MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_FALSE, 1, &nullSpaceVec, &nsp); CHKERRQ(ierr);
  ierr = KSPSetNullSpace(NavierStokesSolver<dim>::ksp2, nsp);
  ierr = MatNullSpaceDestroy(&nsp);

  return 0;
} // setNullSpace