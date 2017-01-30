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

  PetscMPIInt numProcs;
  ierr = MPI_Comm_size(PETSC_COMM_WORLD, &numProcs); CHKERRQ(ierr);

  ierr = setLocalIndexPointsBodies(); CHKERRQ(ierr);
  std::vector<PetscInt> localNumLagrangianPoints(numProcs);
  for (PetscInt procIdx=0; procIdx<numProcs; procIdx++)
  {
    localNumLagrangianPoints[procIdx] = 0;
    for (PetscInt i=0; i<numBodies; i++)
    {
      localNumLagrangianPoints[procIdx] += bodies[i].localNumPoints[procIdx];
    }
  }
  PetscInt globalNumLagrangianPoints = 0;
  for (PetscInt i=0; i<numBodies; i++)
  {
    globalNumLagrangianPoints += bodies[i].numPoints;
  }
  ierr = DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE,
                      globalNumLagrangianPoints, dim, 0, &localNumLagrangianPoints.front(),
                      &bda); CHKERRQ(ierr);
  ierr = PetscObjectViewFromOptions((PetscObject) bda, NULL, "-bda_dmda_view"); CHKERRQ(ierr);
  ierr = DMCompositeAddDM(NavierStokesSolver<dim>::lambdaPack, bda); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // createDMs
