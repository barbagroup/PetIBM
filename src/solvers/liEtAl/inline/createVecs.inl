/*! Implementation of the method `createVecs` of the class `LiEtAlSolver`.
 * \file createVecs.inl
 */


/**
 * \brief Creates the different vectors of the solver.
 */
template <PetscInt dim>
PetscErrorCode LiEtAlSolver<dim>::createVecs()
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  ierr = NavierStokesSolver<dim>::createVecs(); CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(bda, &fTilde); CHKERRQ(ierr);
  ierr = VecDuplicate(fTilde, &rhsf); CHKERRQ(ierr);
  ierr = VecDuplicate(fTilde, &dfTilde); CHKERRQ(ierr);
  ierr = VecDuplicate(NavierStokesSolver<dim>::lambda, &dlambda); CHKERRQ(ierr);
  ierr = VecDuplicate(NavierStokesSolver<dim>::q, &tmp); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // createVecs
