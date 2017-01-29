/*! Implementation of the method `createDMs` of the class `TairaColoniusSolver`.
 * \file createDMs.inl
 */


/*!
 * \brief Creates the DMDA object that accounts for the number of Lagrangian points.
 */
template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::createDMs()
{
  PetscErrorCode ierr;
  
  PetscFunctionBeginUser;

  ierr = NavierStokesSolver<dim>::createDMs(); CHKERRQ(ierr);
  for (PetscInt i=0; i<numBodies; i++)
  {
    ierr = getBodyInfo(bodies[i], numBoundaryPointsOnProcess, boundaryPointIndices); CHKERRQ(ierr);
    ierr = DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, 
                        bodies[i].numPoints, dim, 0, &numBoundaryPointsOnProcess.front(), 
                        &bda); CHKERRQ(ierr);
    ierr = PetscObjectViewFromOptions((PetscObject) bda, NULL, "-bda_dmda_view"); CHKERRQ(ierr);
    ierr = DMCompositeAddDM(NavierStokesSolver<dim>::lambdaPack, bda); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
} // createDMs
