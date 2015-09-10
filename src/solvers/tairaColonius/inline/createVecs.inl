/***************************************************************************//**
 * \file createVecs.inl
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the method `createVecs` of the class `TairaColoniusSolver`.
 */


/**
 * \brief Creates the different vectors of the solver.
 */
template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::createVecs()
{
  PetscErrorCode ierr;

  ierr = NavierStokesSolver<dim>::createVecs();
  ierr = VecDuplicate(NavierStokesSolver<dim>::q, &regularizedForce); CHKERRQ(ierr);
  ierr = VecDuplicate(NavierStokesSolver<dim>::lambda, &nullSpaceVec); CHKERRQ(ierr);

  return 0;
} // createVecs