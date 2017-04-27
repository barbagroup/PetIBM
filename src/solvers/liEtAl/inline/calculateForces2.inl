/*! Methods to calculate the forces acting on a immersed body.
 * \file calculateForces.inl
 */


/*!
 * \brief Computes the forces acting on the immersed body.
 */
template <PetscInt dim>
PetscErrorCode LiEtAlSolver<dim>::calculateForces2()
{
  return 0;
} // calculateForces2


// two-dimensional specialization
template <>
PetscErrorCode LiEtAlSolver<2>::calculateForces2()
{
  PetscErrorCode ierr;

  PetscInt i, j,           // loop indices
           m, n,           // local number of nodes along each direction
           mstart, nstart; // starting indices

  PetscFunctionBeginUser;

  ierr = MatMult(ET, fTilde, tmp); CHKERRQ(ierr);

  Vec fxGlobal, fyGlobal;
  ierr = DMCompositeGetAccess(qPack, tmp, &fxGlobal, &fyGlobal); CHKERRQ(ierr);

  PetscReal localForces[2];

  // compute force in the x-direction
  PetscReal **fx;
  ierr = DMDAVecGetArray(uda, fxGlobal, &fx); CHKERRQ(ierr);
  ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
  localForces[0] = 0.0;
  for (j=nstart; j<nstart+n; j++)
  {
    for (i=mstart; i<mstart+m; i++)
    {
      localForces[0] += mesh->dy[j] * fx[j][i];
    }
  }
  ierr = DMDAVecRestoreArray(uda, fxGlobal, &fx); CHKERRQ(ierr);

  // compute force in the y-direction
  PetscReal **fy;
  ierr = DMDAVecGetArray(vda, fyGlobal, &fy); CHKERRQ(ierr);
  ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
  localForces[1] = 0.0;
  for (j=nstart; j<nstart+n; j++)
  {
    for (i=mstart; i<mstart+m; i++)
    {
      localForces[1] += mesh->dx[i] * fy[j][i];
    }
  }
  ierr = DMDAVecRestoreArray(vda, fyGlobal, &fy); CHKERRQ(ierr);

  ierr = DMCompositeRestoreAccess(qPack, tmp, &fxGlobal, &fyGlobal); CHKERRQ(ierr);

  ierr = MPI_Reduce(localForces, bodyForces, 2, MPIU_REAL, MPI_SUM, 0, MPI_COMM_WORLD); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // calculateForces2


// three-dimensional specialization
template <>
PetscErrorCode LiEtAlSolver<3>::calculateForces2()
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  PetscInt i, j, k,                // loop indices
           m, n, p,                // local number of nodes along each direction
           mstart, nstart, pstart; // starting indices
  
  ierr = MatMult(ET, fTilde, tmp);

  Vec fxGlobal, fyGlobal, fzGlobal;
  ierr = DMCompositeGetAccess(qPack, tmp, &fxGlobal, &fyGlobal, &fzGlobal); CHKERRQ(ierr);
  
  PetscReal localForces[3];

  // compute force in the x-direction
  PetscReal ***fx;
  ierr = DMDAVecGetArray(uda, fxGlobal, &fx); CHKERRQ(ierr);
  ierr = DMDAGetCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  localForces[0] = 0;
  for (k=pstart; k<pstart+p; k++)
  {
    for (j=nstart; j<nstart+n; j++)
    {
      for (i=mstart; i<mstart+m; i++)
      {
        localForces[0] += mesh->dy[j]*mesh->dz[k] * fx[k][j][i];
      }
    }
  }
  ierr = DMDAVecRestoreArray(uda, fxGlobal, &fx); CHKERRQ(ierr);

  // compute force in the y-direction
  PetscReal ***fy;
  ierr = DMDAVecGetArray(vda, fyGlobal, &fy); CHKERRQ(ierr);
  ierr = DMDAGetCorners(vda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  localForces[1] = 0;
  for (k=pstart; k<pstart+p; k++)
  {
    for (j=nstart; j<nstart+n; j++)
    {
      for (i=mstart; i<mstart+m; i++)
      {
        localForces[1] += mesh->dz[k]*mesh->dx[i] * fy[k][j][i];
      }
    }
  }
  ierr = DMDAVecRestoreArray(vda, fyGlobal, &fy); CHKERRQ(ierr);

  // compute force in the z-direction
  PetscReal ***fz;
  ierr = DMDAVecGetArray(wda, fzGlobal, &fz); CHKERRQ(ierr);
  ierr = DMDAGetCorners(wda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  localForces[2] = 0;
  for (k=pstart; k<pstart+p; k++)
  {
    for (j=nstart; j<nstart+n; j++)
    {
      for (i=mstart; i<mstart+m; i++)
      {
        localForces[2] += mesh->dx[i]*mesh->dy[j] * fz[k][j][i];
      }
    }
  }
  ierr = DMDAVecRestoreArray(wda, fzGlobal, &fz); CHKERRQ(ierr);

  ierr = DMCompositeRestoreAccess(qPack, tmp, &fxGlobal, &fyGlobal, &fzGlobal); CHKERRQ(ierr);

  ierr = MPI_Reduce(localForces, bodyForces, 3, MPIU_REAL, MPI_SUM, 0, MPI_COMM_WORLD); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // calculateForces2
