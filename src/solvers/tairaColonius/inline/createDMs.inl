/***************************************************************************//**
 * \file createDMs.inl
 * \author Anush Krishnan (anush@bu.edu), Olivier Mesnard (mesnardo@gwu.edu)
 * \brief Implementation of the method `createDMs` of the class `TairaColoniusSolver`.
 */


/**
 * \brief Creates the DMDA object that accounts for the number of Lagrangian points.
 */
template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::createDMs()
{
  PetscErrorCode ierr;
  
  ierr = NavierStokesSolver<dim>::createDMs(); CHKERRQ(ierr);
  ierr = generateBodyInfo(); CHKERRQ(ierr);
  ierr = DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, 
                      bodies[0].numPoints, dim, 0, &numBoundaryPointsOnProcess.front(), 
                      &bda); CHKERRQ(ierr);
  ierr = DMCompositeAddDM(NavierStokesSolver<dim>::lambdaPack, bda); CHKERRQ(ierr);

  return 0;
} // createDMs