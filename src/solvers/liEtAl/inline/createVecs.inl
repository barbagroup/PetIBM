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
  ierr = VecDuplicate(NavierStokesSolver<dim>::q, &qStar2); CHKERRQ(ierr);
  ierr = VecDuplicate(NavierStokesSolver<dim>::q, &temp); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // createVecs
