/***************************************************************************//**
 * \file generateR2.inl
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the method `generateR2` of the class `TairaColoniusSolver`.
 */


/**
 * \brief Assembles part of the rhs of the Poisson system 
 *        (only values that arise from the boundary ghosts).
 */
template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::generateR2()
{
  return 0;
} // generateR2


// two-dimensional specialization
template <>
PetscErrorCode TairaColoniusSolver<2>::generateR2()
{
  PetscErrorCode ierr;

  PetscInt i, j,           // loop indices 
           m, n,           // local number of nodes along each direction
           mstart, nstart; // starting indices

  ierr = VecSet(r2, 0.0); CHKERRQ(ierr);
  Vec bc2Global;
  ierr = DMCompositeGetAccess(lambdaPack, r2,  &bc2Global, NULL); CHKERRQ(ierr);

  PetscReal **qx, **qy;
  ierr = DMDAVecGetArray(uda, qxLocal, &qx); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(vda, qyLocal, &qy); CHKERRQ(ierr);

  PetscReal **bc2;
  ierr = DMDAVecGetArray(pda, bc2Global, &bc2); CHKERRQ(ierr);
  ierr = DMDAGetCorners(pda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
  // left and right boundaries
  if (flow->boundaries[XPLUS][0].type != PERIODIC) // do not update if x-periodic
  {
    for (j=nstart; j<nstart+n; j++) // loop over y-direction
    {
      if (mstart == 0) // left boundary on current process
      {
        bc2[j][0] -= qx[j][-1];
      }
      if (mstart+m-1 == mesh->nx-1) // right boundary on current process
      {
        bc2[j][mesh->nx-1] += qx[j][mesh->nx-1];
      }
    }
  }
  // bottom and top boundaries
  if (flow->boundaries[YPLUS][1].type != PERIODIC) // do not update if y-periodic
  {
    for (i=mstart; i<mstart+m; i++) // loop over x-direction
    {
      if (nstart == 0) // bottom boundary on current process
      {
        bc2[0][i] -= qy[-1][i];
      }
      if (nstart+n-1 == mesh->ny-1) // top boundary on current process
      {
        bc2[mesh->ny-1][i] += qy[mesh->ny-1][i];
      }
    }
  }
  ierr = DMDAVecRestoreArray(pda, bc2Global, &bc2); CHKERRQ(ierr);

  ierr = DMDAVecRestoreArray(uda, qxLocal, &qx); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(vda, qyLocal, &qy); CHKERRQ(ierr);

  ierr = DMCompositeRestoreAccess(lambdaPack, r2,  &bc2Global, NULL); CHKERRQ(ierr);

  ierr = PetscObjectViewFromOptions((PetscObject) r2, NULL, "-r2_vec_view"); CHKERRQ(ierr);

  return 0;
} // generateR2


// three-dimensional specialization
template <>
PetscErrorCode TairaColoniusSolver<3>::generateR2()
{
  PetscErrorCode ierr;

  PetscInt i, j, k,                // loop indices
           m, n, p,                // local number of nodes along each direction
           mstart, nstart, pstart; // starting indices

  ierr = VecSet(r2, 0.0); CHKERRQ(ierr);
  Vec bc2Global;
  ierr = DMCompositeGetAccess(lambdaPack, r2,  &bc2Global, NULL); CHKERRQ(ierr);

  PetscReal ***qx, ***qy, ***qz;
  ierr = DMDAVecGetArray(uda, qxLocal, &qx); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(vda, qyLocal, &qy); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(wda, qzLocal, &qz); CHKERRQ(ierr);

  PetscReal ***bc2;
  ierr = DMDAVecGetArray(pda, bc2Global, &bc2); CHKERRQ(ierr);
  ierr = DMDAGetCorners(pda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  // left and right boundaries
  if (flow->boundaries[XPLUS][0].type != PERIODIC) // do not update if x-periodic
  {
    for (k=pstart; k<pstart+p; k++) // loop over z-direction
    {
      for (j=nstart; j<nstart+n; j++) // loop over y-direction
      {
        if (mstart == 0) // left boundary on current process
        {
          bc2[k][j][0] -= qx[k][j][-1];
        }
        if (mstart+m-1 == mesh->nx-1) // right boundary on current process
        {
          bc2[k][j][mesh->nx-1] += qx[k][j][mesh->nx-1];
        }
      }
    }
  }
  // bottom and top boundaries
  if (flow->boundaries[YPLUS][1].type != PERIODIC) // do not update if y-periodic
  {
    for (k=pstart; k<pstart+p; k++) // loop over z-direction
    {
      for (i=mstart; i<mstart+m; i++) // loop over x-direction
      {
        if (nstart == 0) // bottom boundary on current process
        {
          bc2[k][0][i] -= qy[k][-1][i];
        }
        if (nstart+n-1 == mesh->ny-1) // top boundary on current process
        {
          bc2[k][mesh->ny-1][i] += qy[k][mesh->ny-1][i];
        }
      }
    }
  }
  // back and front boundaries
  if (flow->boundaries[ZPLUS][2].type != PERIODIC) // do not update if z-periodic
  {
    for (j=nstart; j<nstart+n; j++) // loop over y-direction
    {
      for (i=mstart; i<mstart+m; i++) // loop over x-direction
      {
        if (pstart == 0) // back boundary on current process
        {
          bc2[0][j][i] -= qz[-1][j][i];
        }
        if (pstart+p-1 == mesh->nz-1) // front boundary on current process
        {
          bc2[mesh->nz-1][j][i] += qz[mesh->nz-1][j][i];
        }
      }
    }
  }
  ierr = DMDAVecRestoreArray(pda, bc2Global, &bc2); CHKERRQ(ierr);

  ierr = DMDAVecRestoreArray(uda, qxLocal, &qx); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(vda, qyLocal, &qy); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(wda, qzLocal, &qz); CHKERRQ(ierr);

  ierr = DMCompositeRestoreAccess(lambdaPack, r2,  &bc2Global, NULL); CHKERRQ(ierr);

  ierr = PetscObjectViewFromOptions((PetscObject) r2, NULL, "-r2_vec_view"); CHKERRQ(ierr);

  return 0;
} // generateR2