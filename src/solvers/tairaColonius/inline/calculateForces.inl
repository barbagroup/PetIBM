/*! Methods to calculate the forces acting on a immersed body.
 * \file calculateForces.inl
 */


/*!
 * \brief Computes the forces acting on the immersed bodies.
 * 
 * Sum directly over the Lagrangian forces.
 */
template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::calculateForces()
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  // get access to the global body forces vector
  Vec fGlobal;
  ierr = DMCompositeGetAccess(NavierStokesSolver<dim>::lambdaPack, NavierStokesSolver<dim>::lambda, NULL, &fGlobal); CHKERRQ(ierr);
  PetscReal **f;
  ierr = DMDAVecGetArrayDOF(bda, fGlobal, &f); CHKERRQ(ierr);

  DMDALocalInfo info;
  ierr = DMDAGetLocalInfo(bda, &info); CHKERRQ(ierr);

  PetscReal localForces[dim];
  PetscInt start = 0, end;
  PetscMPIInt rank;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  for (auto &body : bodies)
  {
    end = start + body.numPointsOnProcess[rank];
    for (PetscInt d=0; d<dim; d++)
    {
      localForces[d] = 0.0;
      for (PetscInt i=0; i<info.xm; i++)
      {
        if (start <= i && i < end)
        {
          localForces[d] += f[info.xs + i][d];
        }
      }
    }
    ierr = MPI_Reduce(localForces, body.forces, dim, MPIU_REAL, MPI_SUM, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
    start = end;
  }

  ierr = DMDAVecRestoreArrayDOF(bda, fGlobal, &f); CHKERRQ(ierr);
  ierr = DMCompositeRestoreAccess(NavierStokesSolver<dim>::lambdaPack, NavierStokesSolver<dim>::lambda, NULL, &fGlobal); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // calculateForces
