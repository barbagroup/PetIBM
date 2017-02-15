/*! Implementation of the method `createDMs` of the class `TairaColoniusSolver`.
 * \file createDMs.inl
 */


/*!
 * \brief Creates the DMDA objects for the pressure, fluxes, and Lagrangian forces.
 */
template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::createDMs()
{
  PetscErrorCode ierr;
  
  PetscFunctionBeginUser;

  ierr = NavierStokesSolver<dim>::createDMs(); CHKERRQ(ierr);

  ierr = registerLagPointsOnProcess(); CHKERRQ(ierr);

  // get the total number of Lagrangian points
  PetscInt numLagPoints;
  ierr = getNumLagPoints(numLagPoints); CHKERRQ(ierr);

  // get the number of Lagrangian points per process
  std::vector<PetscInt> numLagPointsOnProcess;
  ierr = getNumLagPointsOnProcess(numLagPointsOnProcess); CHKERRQ(ierr);

  ierr = DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE,
                      numLagPoints, dim, 0, &numLagPointsOnProcess.front(),
                      &bda); CHKERRQ(ierr);
  ierr = PetscObjectViewFromOptions((PetscObject) bda, NULL,
                                    "-bda_dmda_view"); CHKERRQ(ierr);
  ierr = DMCompositeAddDM(NavierStokesSolver<dim>::lambdaPack,
                          bda); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // createDMs
