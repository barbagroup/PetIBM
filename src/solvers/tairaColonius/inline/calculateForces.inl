/*! Methods to calculate the forces acting on a immersed body.
 * \file calculateForces.inl
 */


/**
 * \brief Computes the forces acting on the immersed body.
 */
template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::calculateForces()
{
  return 0;
} // calculateForces


// two-dimensional specialization
template <>
PetscErrorCode TairaColoniusSolver<2>::calculateForces()
{
  PetscErrorCode ierr;

  PetscInt i, j,           // loop indices
           m, n,           // local number of nodes along each direction
           mstart, nstart; // starting indices

  PetscReal forceOnProcess[2];

  // get access to the global body forces vector
  Vec fGlobal;
  ierr = DMCompositeGetAccess(lambdaPack, lambda, NULL, &fGlobal); CHKERRQ(ierr);
  ierr = MatMult(ET, fGlobal, regularizedForce);

  Vec fxGlobal, fyGlobal;
  ierr = DMCompositeGetAccess(qPack, regularizedForce, &fxGlobal, &fyGlobal); CHKERRQ(ierr);

  // force in x-direction
  PetscReal **fx;
  ierr = DMDAVecGetArray(uda, fxGlobal, &fx); CHKERRQ(ierr);
  ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
  forceOnProcess[0] = 0;
  for (j=nstart; j<nstart+n; j++)
  {
    for (i=mstart; i<mstart+m; i++)
    {
      forceOnProcess[0] += mesh->dy[j] * fx[j][i];
    }
  }
  ierr = DMDAVecRestoreArray(uda, fxGlobal, &fx); CHKERRQ(ierr);

  // force in y-direction
  PetscReal **fy;
  ierr = DMDAVecGetArray(vda, fyGlobal, &fy); CHKERRQ(ierr);
  ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
  forceOnProcess[1] = 0;
  for (j=nstart; j<nstart+n; j++)
  {
    for (i=mstart; i<mstart+m; i++)
    {
      forceOnProcess[1] += mesh->dx[i] * fy[j][i];
    }
  }
  ierr = DMDAVecRestoreArray(vda, fyGlobal, &fy); CHKERRQ(ierr);

  ierr = DMCompositeRestoreAccess(qPack, regularizedForce, &fxGlobal, &fyGlobal); CHKERRQ(ierr);
  ierr = DMCompositeRestoreAccess(lambdaPack, lambda, NULL, &fGlobal); CHKERRQ(ierr);

  ierr = MPI_Reduce(forceOnProcess, force, 2, MPIU_REAL, MPI_SUM, 0, MPI_COMM_WORLD); CHKERRQ(ierr);

  return 0;
} // calculateForces


// three-dimensional specialization
template <>
PetscErrorCode TairaColoniusSolver<3>::calculateForces()
{
  PetscErrorCode ierr;

  PetscInt i, j, k,                // loop indices
           m, n, p,                // local number of nodes along each direction
           mstart, nstart, pstart; // starting indices
  
  PetscReal forceOnProcess[3];

  // get access to global body forces vector
  Vec fGlobal;
  ierr = DMCompositeGetAccess(lambdaPack, lambda, NULL, &fGlobal); CHKERRQ(ierr);
  ierr = MatMult(ET, fGlobal, regularizedForce);

  Vec fxGlobal, fyGlobal, fzGlobal;
  ierr = DMCompositeGetAccess(qPack, regularizedForce, &fxGlobal, &fyGlobal, &fzGlobal); CHKERRQ(ierr);
  
  // force in x-direction
  PetscReal ***fx;
  ierr = DMDAVecGetArray(uda, fxGlobal, &fx); CHKERRQ(ierr);
  ierr = DMDAGetCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  forceOnProcess[0] = 0;
  for (k=pstart; k<pstart+p; k++)
  {
    for (j=nstart; j<nstart+n; j++)
    {
      for (i=mstart; i<mstart+m; i++)
      {
        forceOnProcess[0] += mesh->dy[j]*mesh->dz[k] * fx[k][j][i];
      }
    }
  }
  ierr = DMDAVecRestoreArray(uda, fxGlobal, &fx); CHKERRQ(ierr);

  // force in y-direction
  PetscReal ***fy;
  ierr = DMDAVecGetArray(vda, fyGlobal, &fy); CHKERRQ(ierr);
  ierr = DMDAGetCorners(vda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  forceOnProcess[1] = 0;
  for (k=pstart; k<pstart+p; k++)
  {
    for (j=nstart; j<nstart+n; j++)
    {
      for (i=mstart; i<mstart+m; i++)
      {
        forceOnProcess[1] += mesh->dz[k]*mesh->dx[i] * fy[k][j][i];
      }
    }
  }
  ierr = DMDAVecRestoreArray(vda, fyGlobal, &fy); CHKERRQ(ierr);

  // force in z-direction
  PetscReal ***fz;
  ierr = DMDAVecGetArray(wda, fzGlobal, &fz); CHKERRQ(ierr);
  ierr = DMDAGetCorners(wda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  forceOnProcess[2] = 0;
  for (k=pstart; k<pstart+p; k++)
  {
    for (j=nstart; j<nstart+n; j++)
    {
      for (i=mstart; i<mstart+m; i++)
      {
        forceOnProcess[2] += mesh->dx[i]*mesh->dy[j] * fz[k][j][i];
      }
    }
  }
  ierr = DMDAVecRestoreArray(wda, fzGlobal, &fz); CHKERRQ(ierr);

  ierr = DMCompositeRestoreAccess(qPack, regularizedForce, &fxGlobal, &fyGlobal, &fzGlobal); CHKERRQ(ierr);
  ierr = DMCompositeRestoreAccess(lambdaPack, lambda, NULL, &fGlobal); CHKERRQ(ierr);

  ierr = MPI_Reduce(forceOnProcess, force, 3, MPIU_REAL, MPI_SUM, 0, MPI_COMM_WORLD); CHKERRQ(ierr);

  return 0;
} // calculateForces


/**
 * \brief Computes the forces acting on the immersed body.
 * 
 * Sum directly over the Lagrangian forces.
 */
template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::calculateForcesTC()
{
  PetscErrorCode ierr;

  // get access to the global body forces vector
  Vec fGlobal;
  ierr = DMCompositeGetAccess(NavierStokesSolver<dim>::lambdaPack, NavierStokesSolver<dim>::lambda, NULL, &fGlobal); CHKERRQ(ierr);
  PetscReal **f;
  ierr = DMDAVecGetArrayDOF(bda, fGlobal, &f); CHKERRQ(ierr);
  PetscInt m, mstart;  // local number of points and starting index
  ierr = DMDAGetCorners(bda, &mstart, NULL, NULL, &m, NULL, NULL); CHKERRQ(ierr);
  
  // sum force components over the local Lagrangian points
  PetscReal forceOnProcess[dim];  // local force vector
  for (PetscInt d=0; d<dim; d++)
  {
    forceOnProcess[d] = 0.0;
    for (PetscInt i=mstart; i<mstart+m; i++)
    {
      forceOnProcess[d] += f[i][d];
    }
  }

  ierr = DMDAVecRestoreArrayDOF(bda, fGlobal, &f); CHKERRQ(ierr);
  ierr = DMCompositeRestoreAccess(NavierStokesSolver<dim>::lambdaPack, NavierStokesSolver<dim>::lambda, NULL, &fGlobal); CHKERRQ(ierr);

  ierr = MPI_Reduce(forceOnProcess, force, dim, MPIU_REAL, MPI_SUM, 0, MPI_COMM_WORLD); CHKERRQ(ierr);

  return 0;
} // calculateForcesTC
