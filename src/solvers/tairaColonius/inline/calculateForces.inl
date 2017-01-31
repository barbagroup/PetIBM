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
  PetscInt m, mstart;  // local number of points and starting index
  ierr = DMDAGetCorners(bda, &mstart, NULL, NULL, &m, NULL, NULL); CHKERRQ(ierr);
  
  // sum force components over the local Lagrangian points
  PetscReal localForces[dim];  // local force vector
  PetscInt start = 0, end;
  for (PetscInt bIdx=0; bIdx<numBodies; bIdx++)
  {
    end = start + bodies[bIdx].numPoints;
    for (PetscInt d=0; d<dim; d++)
    {
      localForces[d] = 0.0;
      for (PetscInt i=mstart; i<mstart+m; i++)
      {
        if (start <= i && i < end)
        {
          localForces[d] += f[i][d];
        }
      }
    }
    ierr = MPI_Reduce(localForces, bodies[bIdx].forces, dim, MPIU_REAL, MPI_SUM, 0, MPI_COMM_WORLD); CHKERRQ(ierr);
    start = end;
  }

  ierr = DMDAVecRestoreArrayDOF(bda, fGlobal, &f); CHKERRQ(ierr);
  ierr = DMCompositeRestoreAccess(NavierStokesSolver<dim>::lambdaPack, NavierStokesSolver<dim>::lambda, NULL, &fGlobal); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // calculateForces
