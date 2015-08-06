/***************************************************************************//**
 * \file initializeLambda.inl
 * \author Anush Kirshnan (anush@bu.edu), Olivier Mesnard (mesnardo@gwu.edu)
 * \brief Implementation of the method to initialize the pressure field.
 */

/**
 * \brief Initiliazes the lambda vector that contains the pressure.
 *
 * Only reads the pressure field when the simulation is restarted or when a 
 * custom initial field is used.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::initializeLambda()
{
  PetscErrorCode ierr;

  if (simParams->startStep > 0 || flowDesc->initialCustomField)
  {
    ierr = readLambda(); CHKERRQ(ierr);
  }

  return 0;
} // initializeLambda