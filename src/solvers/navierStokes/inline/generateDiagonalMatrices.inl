/***************************************************************************//**
 * \file generateDiagonalMatrices.inl
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the method `generateDiagonalMatrices` 
 *        of the class `NavierStokesSolver`.
 */


/**
 * \brief Assembles the diagonal matrices \f$ B^N \f$, \f$ \hat{M} \f$ 
 *        and \f$ R^{-1} \f$.
 *
 * The matrix \f$ B^N \f$ is the approximate inverse of the the matrix \f$ A \f$ 
 * resulting from the implicit contributions in the momentum equations.
 * The first-order approximation is a diagonal matrix.
 * The matrices \f$ \hat{M} \f$ and \f$ R^{-1} \f$ are used to scale 
 * the Stokes system to make it symmetric.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::generateDiagonalMatrices()
{
  return 0;
} // generateDiagonalMatrices


// two-dimensional specialization
template <>
PetscErrorCode NavierStokesSolver<2>::generateDiagonalMatrices()
{
  PetscErrorCode ierr;

  PetscInt i, j,           // loop indices
           m, n,           // local number of nodes along each direction
           mstart, nstart; // starting indices

  Vec MHatxGlobal, MHatyGlobal;
  ierr = DMCompositeGetAccess(qPack, MHat, &MHatxGlobal, &MHatyGlobal); CHKERRQ(ierr);
  Vec RInvxGlobal, RInvyGlobal;
  ierr = DMCompositeGetAccess(qPack, RInv, &RInvxGlobal, &RInvyGlobal); CHKERRQ(ierr);
  Vec BNxGlobal, BNyGlobal;
  ierr = DMCompositeGetAccess(qPack, BN, &BNxGlobal, &BNyGlobal); CHKERRQ(ierr);

  // x-direction
  PetscReal **MHatx;
  ierr = DMDAVecGetArray(uda, MHatxGlobal, &MHatx); CHKERRQ(ierr);
  PetscReal **RInvx;
  ierr = DMDAVecGetArray(uda, RInvxGlobal, &RInvx); CHKERRQ(ierr);
  PetscReal **BNx;
  ierr = DMDAVecGetArray(uda, BNxGlobal, &BNx); CHKERRQ(ierr);
  ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
  for (j=nstart; j<nstart+n; j++)
  {
    for (i=mstart; i<mstart+m; i++)
    {
      MHatx[j][i] = (i<mesh->nx-1) ? 0.5*(mesh->dx[i]+mesh->dx[i+1]) : 0.5*(mesh->dx[i]+mesh->dx[0]);
      RInvx[j][i] = 1.0/mesh->dy[j];
      BNx[j][i]   = parameters->dt/(MHatx[j][i]*RInvx[j][i]);
    }
  }
  ierr = DMDAVecRestoreArray(uda, MHatxGlobal, &MHatx); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(uda, RInvxGlobal, &RInvx); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(uda, BNxGlobal, &BNx); CHKERRQ(ierr);

  // y-direction
  PetscReal **MHaty;
  ierr = DMDAVecGetArray(vda, MHatyGlobal, &MHaty); CHKERRQ(ierr);
  PetscReal **RInvy;
  ierr = DMDAVecGetArray(vda, RInvyGlobal, &RInvy); CHKERRQ(ierr);
  PetscReal **BNy;
  ierr = DMDAVecGetArray(vda, BNyGlobal, &BNy); CHKERRQ(ierr);
  ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
  for (j=nstart; j<nstart+n; j++)
  {
    for (i=mstart; i<mstart+m; i++)
    {
      MHaty[j][i] = (j<mesh->ny-1) ? 0.5*(mesh->dy[j]+mesh->dy[j+1]) : 0.5*(mesh->dy[j]+mesh->dy[0]);
      RInvy[j][i] = 1.0/mesh->dx[i];
      BNy[j][i] = parameters->dt/(MHaty[j][i]*RInvy[j][i]);
    }
  }
  ierr = DMDAVecRestoreArray(vda, MHatyGlobal, &MHaty); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(vda, RInvyGlobal, &RInvy); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(vda, BNyGlobal, &BNy); CHKERRQ(ierr);

  ierr = DMCompositeRestoreAccess(qPack, MHat, &MHatxGlobal, &MHatyGlobal); CHKERRQ(ierr);
  ierr = DMCompositeRestoreAccess(qPack, RInv, &RInvxGlobal, &RInvyGlobal); CHKERRQ(ierr);
  ierr = DMCompositeRestoreAccess(qPack, BN, &BNxGlobal, &BNyGlobal); CHKERRQ(ierr);

  return 0;
} // generateDiagonalMatrices


// three-dimensional specialization
template <>
PetscErrorCode NavierStokesSolver<3>::generateDiagonalMatrices()
{
  PetscErrorCode ierr;

  PetscInt i, j, k,                // loop indices
           m, n, p,                // local number of nodes along each direction
           mstart, nstart, pstart; // starting indices
  
  Vec MHatxGlobal, MHatyGlobal, MHatzGlobal;
  ierr = DMCompositeGetAccess(qPack, MHat, &MHatxGlobal, &MHatyGlobal, &MHatzGlobal); CHKERRQ(ierr);
  Vec RInvxGlobal, RInvyGlobal, RInvzGlobal;
  ierr = DMCompositeGetAccess(qPack, RInv, &RInvxGlobal, &RInvyGlobal, &RInvzGlobal); CHKERRQ(ierr);
  Vec BNxGlobal, BNyGlobal, BNzGlobal;
  ierr = DMCompositeGetAccess(qPack, BN, &BNxGlobal, &BNyGlobal, &BNzGlobal); CHKERRQ(ierr);
  
  // x-direction
  PetscReal ***MHatx;
  ierr = DMDAVecGetArray(uda, MHatxGlobal, &MHatx); CHKERRQ(ierr);
  PetscReal ***RInvx;
  ierr = DMDAVecGetArray(uda, RInvxGlobal, &RInvx); CHKERRQ(ierr);
  PetscReal ***BNx;
  ierr = DMDAVecGetArray(uda, BNxGlobal, &BNx); CHKERRQ(ierr);
  ierr = DMDAGetCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  for (k=pstart; k<pstart+p; k++)
  {
    for (j=nstart; j<nstart+n; j++)
    {
      for (i=mstart; i<mstart+m; i++)
      {
        MHatx[k][j][i] = (i < mesh->nx-1) ? 0.5*(mesh->dx[i]+mesh->dx[i+1]) : 0.5*(mesh->dx[i]+mesh->dx[0]);
        RInvx[k][j][i] = 1.0/(mesh->dy[j]*mesh->dz[k]);
        BNx[k][j][i]   = parameters->dt/(MHatx[k][j][i]*RInvx[k][j][i]);
      }
    }
  }
  ierr = DMDAVecRestoreArray(uda, MHatxGlobal, &MHatx); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(uda, RInvxGlobal, &RInvx); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(uda, BNxGlobal, &BNx); CHKERRQ(ierr);
  
  // y-direction
  PetscReal ***MHaty;
  ierr = DMDAVecGetArray(vda, MHatyGlobal, &MHaty); CHKERRQ(ierr);
  PetscReal ***RInvy;
  ierr = DMDAVecGetArray(vda, RInvyGlobal, &RInvy); CHKERRQ(ierr);
  PetscReal ***BNy;
  ierr = DMDAVecGetArray(vda, BNyGlobal, &BNy); CHKERRQ(ierr);
  ierr = DMDAGetCorners(vda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  for (k=pstart; k<pstart+p; k++)
  {
    for (j=nstart; j<nstart+n; j++)
    {
      for (i=mstart; i<mstart+m; i++)
      {
        MHaty[k][j][i] = (j < mesh->ny-1) ? 0.5*(mesh->dy[j]+mesh->dy[j+1]) : 0.5*(mesh->dy[j]+mesh->dy[0]);
        RInvy[k][j][i] = 1.0/(mesh->dx[i]*mesh->dz[k]);
        BNy[k][j][i]   = parameters->dt/(MHaty[k][j][i]*RInvy[k][j][i]);
      }
    }
  }
  ierr = DMDAVecRestoreArray(vda, MHatyGlobal, &MHaty); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(vda, RInvyGlobal, &RInvy); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(vda, BNyGlobal, &BNy); CHKERRQ(ierr);
  
  // z-direction
  PetscReal ***MHatz;
  ierr = DMDAVecGetArray(wda, MHatzGlobal, &MHatz); CHKERRQ(ierr);
  PetscReal ***RInvz;
  ierr = DMDAVecGetArray(wda, RInvzGlobal, &RInvz); CHKERRQ(ierr);
  PetscReal ***BNz;
  ierr = DMDAVecGetArray(wda, BNzGlobal, &BNz); CHKERRQ(ierr);
  ierr = DMDAGetCorners(wda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  for (k=pstart; k<pstart+p; k++)
  {
    for (j=nstart; j<nstart+n; j++)
    {
      for (i=mstart; i<mstart+m; i++)
      {
        MHatz[k][j][i] = (k < mesh->nz-1) ? 0.5*(mesh->dz[k]+mesh->dz[k+1]) : 0.5*(mesh->dz[k]+mesh->dz[0]);
        RInvz[k][j][i] = 1.0/(mesh->dx[i]*mesh->dy[j]);
        BNz[k][j][i]   = parameters->dt/(MHatz[k][j][i]*RInvz[k][j][i]);
      }
    }
  }
  ierr = DMDAVecRestoreArray(wda, MHatzGlobal, &MHatz); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(wda, RInvzGlobal, &RInvz); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(wda, BNzGlobal, &BNz); CHKERRQ(ierr);
  
  ierr = DMCompositeRestoreAccess(qPack, MHat, &MHatxGlobal, &MHatyGlobal, &MHatzGlobal); CHKERRQ(ierr);
  ierr = DMCompositeRestoreAccess(qPack, RInv, &RInvxGlobal, &RInvyGlobal, &RInvzGlobal); CHKERRQ(ierr);
  ierr = DMCompositeRestoreAccess(qPack, BN, &BNxGlobal, &BNyGlobal, &BNzGlobal); CHKERRQ(ierr);

  return 0;
} // generateDiagonalMatrices