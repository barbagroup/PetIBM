/***************************************************************************//**
 * \file initializeLambda.inl
 * \author Anush Kirshnan (anush@bu.edu), Olivier Mesnard (mesnardo@gwu.edu)
 * \brief Implementation of the method `initializeLambda` 
 *        of the class `NavierStokesSolver`.
 */

/**
 * \brief Initializes the lambda vector that contains the pressure with zeros.
 *
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::initializeLambda()
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  ierr = VecSet(lambda, 0.0); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // initializeLambda
