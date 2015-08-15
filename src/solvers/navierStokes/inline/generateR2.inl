/***************************************************************************//**
 * \file generateR2.inl
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the method `generateR2` of the class `NavierStokesSolver`.
 */


/**
 * \brief
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::generateR2()
{
  return 0;
} // generateR2


// two-dimensional specialization
template <>
PetscErrorCode NavierStokesSolver<2>::generateR2()
{
  PetscErrorCode ierr;

  PetscInt i, j,           // loop indices
           m, n,           // local number of nodes along each direction
           mstart, nstart; // starting indices

  ierr = VecSet(r2, 0.0); CHKERRQ(ierr);
  Vec bc2Global;
  ierr = DMCompositeGetAccess(lambdaPack, r2,  &bc2Global); CHKERRQ(ierr);

  PetscReal **qx, **qy;
  ierr = DMDAVecGetArray(uda, qxLocal, &qx); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(vda, qyLocal, &qy); CHKERRQ(ierr);

  PetscReal **bc2;
  ierr = DMDAVecGetArray(pda, bc2Global, &bc2); CHKERRQ(ierr);
  ierr = DMDAGetCorners(pda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
  // x-faces
  if (flow->boundaries[XPLUS][0].type != PERIODIC) // don't update if the BC type is periodic
  {
    //loop over all points on the x-face
    for (j=nstart; j<nstart+n; j++)
    {
      // -X
      if (mstart == 0) // if the -X face is in the current process
      {
        bc2[j][0] -= qx[j][-1];
      }
      // +X
      if (mstart+m-1 == mesh->nx-1) // if the +X face is in the current process
      {
        bc2[j][mesh->nx-1] += qx[j][mesh->nx-1];
      }
    }
  }
  // y-faces
  if (flow->boundaries[YPLUS][1].type != PERIODIC) // don't update if the BC type is periodic
  {
    //loop over all points on the y-face
    for (i=mstart; i<mstart+m; i++)
    {
      // -Y
      if (nstart == 0) // if the -Y face is in the current process
      {
        bc2[0][i] -= qy[-1][i];
      }
      // +Y
      if (nstart+n-1 == mesh->ny-1) // if the +Y face is in the current process
      {
        bc2[mesh->ny-1][i] += qy[mesh->ny-1][i];
      }
    }
  }
  ierr = DMDAVecRestoreArray(pda, bc2Global, &bc2); CHKERRQ(ierr);

  ierr = DMDAVecRestoreArray(uda, qxLocal, &qx); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(vda, qyLocal, &qy); CHKERRQ(ierr);

  ierr = DMCompositeRestoreAccess(lambdaPack, r2,  &bc2Global); CHKERRQ(ierr);

  return 0;
} // generateR2


// three-dimensional specialization
template <>
PetscErrorCode NavierStokesSolver<3>::generateR2()
{
  PetscErrorCode ierr;
  PetscInt i, j, k, 
           m, n, p, 
           mstart, nstart, pstart;

  ierr = VecSet(r2, 0.0); CHKERRQ(ierr);
  Vec bc2Global;
  ierr = DMCompositeGetAccess(lambdaPack, r2,  &bc2Global); CHKERRQ(ierr);

  PetscReal ***qx, ***qy, ***qz;
  ierr = DMDAVecGetArray(uda, qxLocal, &qx); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(vda, qyLocal, &qy); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(wda, qzLocal, &qz); CHKERRQ(ierr);

  PetscReal ***bc2;
  ierr = DMDAVecGetArray(pda, bc2Global, &bc2); CHKERRQ(ierr);
  ierr = DMDAGetCorners(pda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  // x-faces
  if (flow->boundaries[XPLUS][0].type != PERIODIC) // don't update if the BC type is periodic
  {
    //loop over all points on the x-face
    for (k=pstart; k<pstart+p; k++)
    {
      for (j=nstart; j<nstart+n; j++)
      {
        // -X
        if (mstart == 0) // if the -X face is in the current process
        {
          bc2[k][j][0] -= qx[k][j][-1];
        }
        // +X
        if (mstart+m-1 == mesh->nx-1) // if the +X face is in the current process
        {
          bc2[k][j][mesh->nx-1] += qx[k][j][mesh->nx-1];
        }
      }
    }
  }
  // y-faces
  if (flow->boundaries[YPLUS][1].type != PERIODIC) // don't update if the BC type is periodic
  {
    //loop over all points on the y-face
    for (k=pstart; k<pstart+p; k++)
    {
      for (i=mstart; i<mstart+m; i++)
      {
        // -Y
        if (nstart == 0) // if the -Y face is in the current process
        {
          bc2[k][0][i] -= qy[k][-1][i];
        }
        // +Y
        if (nstart+n-1 == mesh->ny-1) // if the +Y face is in the current process
        {
          bc2[k][mesh->ny-1][i] += qy[k][mesh->ny-1][i];
        }
      }
    }
  }
  // z-faces
  if (flow->boundaries[ZPLUS][2].type != PERIODIC) // don't update if the BC type is periodic
  {
    //loop over all points on the z-face
    for (j=nstart; j<nstart+n; j++)
    {
      for (i=mstart; i<mstart+m; i++)
      {
        // -Z
        if (pstart == 0) // if the -Z face is in the current process
        {
          bc2[0][j][i] -= qz[-1][j][i];
        }
        // +Z
        if (pstart+p-1 == mesh->nz-1) // if the +Z face is in the current process
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

  ierr = DMCompositeRestoreAccess(lambdaPack, r2,  &bc2Global); CHKERRQ(ierr);

  return 0;
} // generateR2