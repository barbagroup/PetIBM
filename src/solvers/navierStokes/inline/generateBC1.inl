/***************************************************************************//**
 * \file generateBC1.inl
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the method `generateBC1` of the class `NavierStokesSolver`.
 */


/**
 * \brief
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::generateBC1()
{
  return 0;
} // generateBC1


// two-dimensional specialization
template <>
PetscErrorCode NavierStokesSolver<2>::generateBC1()
{
  PetscErrorCode ierr;
  PetscInt i, j,           // loop indices
           M, N,           // global number of nodes along each direction
           m, n,           // local number of nodes along each direction
           mstart, nstart; // starting indices
  
  PetscReal nu = flow->nu; // viscosity
  PetscReal alpha = parameters->diffusion.coefficients[0]; // implicit diffusion coefficient
  PetscReal coeffMinus = 0.0, coeffPlus = 0.0;

  PetscReal *dx = &mesh->dx[0],
            *dy = &mesh->dy[0];
  PetscInt nx = mesh->nx,
           ny = mesh->ny;

  ierr = VecSet(bc1, 0.0); CHKERRQ(ierr);
  Vec bc1xGlobal, bc1yGlobal;
  ierr = DMCompositeGetAccess(qPack, bc1, &bc1xGlobal, &bc1yGlobal); CHKERRQ(ierr);
                 
  // fluxes in x-direction
  PetscReal **bc1x;
  ierr = DMDAVecGetArray(uda, bc1xGlobal, &bc1x); CHKERRQ(ierr);
  PetscReal **qx;
  ierr = DMDAVecGetArray(uda, qxLocal, &qx); CHKERRQ(ierr);
  ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
  ierr = DMDAGetInfo(uda, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
  // left and right boundaries
  if (flow->boundaries[XPLUS][0].type != PERIODIC) // do not update if x-periodic
  {
    dxMinus = dx[0];
    dxPlus = dx[1];
    coeffMinus = alpha*nu * 2.0/dxMinus/(dxMinus+dxPlus);
    dxMinus = dx[M-1];
    dxPlus dx
    coeffPlus  = alpha*nu * 2.0/dx[M]/(dx[M-1]+dx[M]);
    for (j=nstart; j<nstart+n; j++) // loop over y-direction
    {
      if (mstart == 0) // left boundary on current process
      {
        switch (flow->boundaries[XMINUS][0].type)
        {
          case DIRICHLET:
            bc1x[j][0] += coeffMinus*qx[j][-1]/dy[j];
            break;
          case CONVECTIVE:
          default:
            break;
        }
      }
      if (mstart+m-1 == M-1) // right on current process
      {     
        switch (flow->boundaries[XPLUS][0].type)
        {
          case DIRICHLET:
            bc1x[j][M-1] += coeffPlus*qx[j][M]/dy[j];
            break;
          case CONVECTIVE:
          default:
            break;
        }
      }
    }
  }
  // bottom and top boundaries
  if (flow->boundaries[YPLUS][0].type != PERIODIC) // do not update if y-periodic
  {
    coeffMinus = alpha*nu*2.0/dyU[0]/(dyU[0]+dyU[1]);
    coeffPlus  = alpha*nu*2.0/dyU[N]/(dyU[N]+dyU[N-1]);
    // loop over all points on the y-face
    for (i=mstart; i<mstart+m; i++)
    { 
      // -Y
      if (nstart == 0) // if the -Y boundary is in the current process
      {
        switch (flow->boundaries[YMINUS][0].type)
        {
          case CONVECTIVE:
          case DIRICHLET : bc1x[0][i] += coeffMinus*qx[-1][i]; break;
          default        : break;
        }
      }
      // +Y
      if (nstart+n-1 == N-1) // if the +Y boundary is in the current process
      {
        switch (flow->boundaries[YPLUS][0].type)
        {
          case CONVECTIVE:
          case DIRICHLET : bc1x[N-1][i] += coeffPlus*qx[N][i]; break;
          default        : break;
        }
      }
    }
  }
  ierr = DMDAVecRestoreArray(uda, bc1xGlobal, &bc1x); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(uda, qxLocal, &qx); CHKERRQ(ierr);
  
  // V-FLUXES
  PetscReal **bc1y;
  ierr = DMDAVecGetArray(vda, bc1yGlobal, &bc1y); CHKERRQ(ierr);
  PetscReal **qy;
  ierr = DMDAVecGetArray(vda, qyLocal, &qy); CHKERRQ(ierr);
  ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
  ierr = DMDAGetInfo(vda, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
  // x-faces
  if (flow->boundaries[XPLUS][1].type != PERIODIC) // don't update if the BC type is periodic
  {
    coeffMinus = alphaImplicit*nu*2.0/dxV[0]/(dxV[0]+dxV[1]);
    coeffPlus  = alphaImplicit*nu*2.0/dxV[M]/(dxV[M]+dxV[M-1]);
    // loop over all points on the x-face
    for (j=nstart; j<nstart+n; j++)
    {
      // -X
      if (mstart == 0)
      {
        switch (flow->boundaries[XMINUS][1].type)
        {
          case CONVECTIVE:
          case DIRICHLET : bc1y[j][0] += coeffMinus*qy[j][-1]; break;
          default        : break;
        }
      }
      // +X
      if (mstart+m-1 == M-1)
      {
        switch (flow->boundaries[XPLUS][1].type)
        {
          case CONVECTIVE:
          case DIRICHLET : bc1y[j][M-1] += coeffPlus*qy[j][M]; break;
          default        : break;
        }
      }
    }
  }
  // y-faces
  if (flow->boundaries[YPLUS][1].type != PERIODIC) // don't update if the BC type is periodic
  {
    coeffMinus = alphaImplicit*nu*2.0/dyV[0]/(dyV[0]+dyV[1]);
    coeffPlus  = alphaImplicit*nu*2.0/dyV[N]/(dyV[N]+dyV[N-1]);
    // loop over all points on the y-face
    for (i=mstart; i<mstart+m; i++)
    { 
      // -Y
      if (nstart == 0)
      {
        switch (flow->boundaries[YMINUS][1].type)
        {
          case CONVECTIVE:
          case DIRICHLET : bc1y[0][i] += coeffMinus*qy[-1][i]/mesh->dx[i]; break;
          default        : break;
        }
      }
      // +Y
      if (nstart+n-1 == N-1)
      {
        switch (flow->boundaries[YPLUS][1].type)
        {
          case CONVECTIVE:
          case DIRICHLET : bc1y[N-1][i] += coeffPlus*qy[N][i]/mesh->dx[i]; break;
          default        : break;
        }
      }
    }
  }
  ierr = DMDAVecRestoreArray(vda, bc1yGlobal, &bc1y); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(vda, qyLocal, &qy); CHKERRQ(ierr);

  ierr = DMCompositeRestoreAccess(qPack, bc1,  &bc1xGlobal, &bc1yGlobal); CHKERRQ(ierr);

  return 0;
} // generateBC1


// three-dimensional specialization
template <>
PetscErrorCode NavierStokesSolver<3>::generateBC1()
{
  PetscErrorCode ierr;
  PetscInt i, j, k, 
           M, N, P, 
           m, n, p, 
           mstart, nstart, pstart;
  
  PetscReal nu = flow->nu;
  PetscReal alphaImplicit = parameters->diffusion.coefficients[0];
  PetscReal coeffMinus = 0.0, coeffPlus = 0.0;

  ierr = VecSet(bc1, 0.0); CHKERRQ(ierr);
  Vec bc1xGlobal, bc1yGlobal, bc1zGlobal;
  ierr = DMCompositeGetAccess(qPack, bc1, &bc1xGlobal, &bc1yGlobal, &bc1zGlobal); CHKERRQ(ierr);
  
  // U-FLUXES
  PetscReal ***bc1x;
  ierr = DMDAVecGetArray(uda, bc1xGlobal, &bc1x); CHKERRQ(ierr);
  PetscReal ***qx;
  ierr = DMDAVecGetArray(uda, qxLocal, &qx); CHKERRQ(ierr);
  ierr = DMDAGetCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  ierr = DMDAGetInfo(uda, NULL, &M, &N, &P, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
  // x-faces
  if (flow->boundaries[XPLUS][0].type != PERIODIC) // don't update if the BC type is periodic
  {
    coeffMinus = alphaImplicit*nu*2.0/dxU[0]/(dxU[0]+dxU[1]);
    coeffPlus  = alphaImplicit*nu*2.0/dxU[M]/(dxU[M]+dxU[M-1]);
    //loop over all points on the x-face
    for (k=pstart; k<pstart+p; k++)
    {
      for (j=nstart; j<nstart+n; j++)
      {
        // -X
        if (mstart == 0) // if the -X face is in the current process
        {
          switch (flow->boundaries[XMINUS][0].type)
          {
            case DIRICHLET :
            case CONVECTIVE: bc1x[k][j][0] += coeffMinus*qx[k][j][-1]/(mesh->dy[j]*mesh->dz[k]); break;
            default        : break;
          }
        }
        // +X
        if (mstart+m-1 == M-1) // if the +X face is in the current process
        {
          switch (flow->boundaries[XPLUS][0].type)
          {
            case DIRICHLET :
            case CONVECTIVE: bc1x[k][j][M-1] += coeffPlus*qx[k][j][M]/(mesh->dy[j]*mesh->dz[k]); break;
            default        : break;
          }
        }
      }
    }
  }
  // y-faces
  if (flow->boundaries[YPLUS][0].type != PERIODIC) // don't update if the BC type is periodic
  {
    coeffMinus = alphaImplicit*nu*2.0/dyU[0]/(dyU[0]+dyU[1]);
    coeffPlus  = alphaImplicit*nu*2.0/dyU[N]/(dyU[N]+dyU[N-1]);
    // loop over all points on the y-face
    for (k=pstart; k<pstart+p; k++)
    {
      for (i=mstart; i<mstart+m; i++)
      {
        // -Y
        if (nstart == 0) // if the -Y boundary is in the current process
        {
          switch (flow->boundaries[YMINUS][0].type)
          {
            case DIRICHLET :
            case CONVECTIVE: bc1x[k][0][i] += coeffMinus*qx[k][-1][i]; break;
            default        : break;
          }
        }
        // +Y
        if (nstart+n-1 == N-1) // if the +Y boundary is in the current process
        {
          switch (flow->boundaries[YPLUS][0].type)
          {
            case DIRICHLET :
            case CONVECTIVE: bc1x[k][N-1][i] += coeffPlus*qx[k][N][i]; break;
            default        : break;
          }
        }
      }
    }
  }
  // z-faces
  if (flow->boundaries[ZPLUS][0].type != PERIODIC) // don't update if the BC type is periodic
  {
    coeffMinus = alphaImplicit*nu*2.0/dzU[0]/(dzU[0]+dzU[1]);
    coeffPlus  = alphaImplicit*nu*2.0/dzU[P]/(dzU[P]+dzU[P-1]);
    // loop over all points on the z-face
    for (j=nstart; j<nstart+n; j++)
    {
      for (i=mstart; i<mstart+m; i++)
      {
        // -Z
        if (pstart == 0) // if the -Z boundary is in the current process
        {
          switch (flow->boundaries[ZMINUS][0].type)
          {
            case DIRICHLET :
            case CONVECTIVE: bc1x[0][j][i] += coeffMinus*qx[-1][j][i]; break;
            default        : break;
          }
        }
        // +Z
        if (pstart+p-1 == P-1) // if the +Z boundary is in the current process
        {
          switch (flow->boundaries[ZPLUS][0].type)
          {
            case DIRICHLET :
            case CONVECTIVE: bc1x[P-1][j][i] += coeffPlus*qx[P][j][i]; break;
            default        : break;
          }
        }
      }
    }
  }
  ierr = DMDAVecRestoreArray(uda, bc1xGlobal, &bc1x); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(uda, qxLocal, &qx); CHKERRQ(ierr);
  
  // V-FLUXES
  PetscReal ***bc1y;
  ierr = DMDAVecGetArray(vda, bc1yGlobal, &bc1y); CHKERRQ(ierr);
  PetscReal ***qy;
  ierr = DMDAVecGetArray(vda, qyLocal, &qy); CHKERRQ(ierr);
  ierr = DMDAGetCorners(vda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  ierr = DMDAGetInfo(vda, NULL, &M, &N, &P, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
  // x-faces
  if (flow->boundaries[XPLUS][1].type != PERIODIC) // don't update if the BC type is periodic
  {
    coeffMinus = alphaImplicit*nu*2.0/dxV[0]/(dxV[0]+dxV[1]);
    coeffPlus  = alphaImplicit*nu*2.0/dxV[M]/(dxV[M]+dxV[M-1]);
    // loop over all points on the x-face
    for (k=pstart; k<pstart+p; k++)
    {
      for (j=nstart; j<nstart+n; j++)
      {
        // -X
        if (mstart == 0) // if the -X face is in the current process
        {
          switch (flow->boundaries[XMINUS][1].type)
          {
            case DIRICHLET :
            case CONVECTIVE: bc1y[k][j][0] += coeffMinus*qy[k][j][-1]; break;
            default        : break;
          }
        }
        // +X
        if (mstart+m-1 == M-1) // if the +X face is in the current process
        {
          switch (flow->boundaries[XPLUS][1].type)
          {
            case DIRICHLET :
            case CONVECTIVE: bc1y[k][j][M-1] += coeffPlus*qy[k][j][M]; break;
            default        : break;
          }
        }
      }
    }
  }
  // y-faces
  if (flow->boundaries[YPLUS][1].type != PERIODIC) // don't update if the BC type is periodic
  {
    coeffMinus = alphaImplicit*nu*2.0/dyV[0]/(dyV[0]+dyV[1]);
    coeffPlus  = alphaImplicit*nu*2.0/dyV[N]/(dyV[N]+dyV[N-1]);
    // loop over all points on the y-face
    for (k=pstart; k<pstart+p; k++)
    {
      for (i=mstart; i<mstart+m; i++)
      {
        // -Y
        if (nstart == 0) // if the -Y face is in the current process
        {
          switch (flow->boundaries[YMINUS][1].type)
          {
            case DIRICHLET :
            case CONVECTIVE: bc1y[k][0][i] += coeffMinus*qy[k][-1][i]/(mesh->dz[k]*mesh->dx[i]); break;
            default        : break;
          }
        }
        // +Y
        if (nstart+n-1 == N-1) // if the +Y face is in the current process
        {
          switch (flow->boundaries[YPLUS][1].type)
          {
            case DIRICHLET :
            case CONVECTIVE: bc1y[k][N-1][i] += coeffPlus*qy[k][N][i]/(mesh->dz[k]*mesh->dx[i]); break;
            default        : break;
          }
        }
      }
    }
  }
  // z-faces
  if (flow->boundaries[ZPLUS][1].type != PERIODIC) // don't update if the BC type is periodic
  {
    coeffMinus = alphaImplicit*nu*2.0/dzV[0]/(dzV[0]+dzV[1]);
    coeffPlus  = alphaImplicit*nu*2.0/dzV[P]/(dzV[P]+dzV[P-1]);
    // loop over all points on the z-face
    for (j=nstart; j<nstart+n; j++)
    {
      for (i=mstart; i<mstart+m; i++)
      {
        // -Z
        if (pstart == 0) // if the -Z face is in the current process
        {
          switch (flow->boundaries[ZMINUS][1].type)
          {
            case DIRICHLET :
            case CONVECTIVE: bc1y[0][j][i] += coeffMinus*qy[-1][j][i]; break;
            default        : break;
          }
        }
        // +Z
        if (pstart+p-1 == P-1) // if the +Z face is in the current process
        {
          switch (flow->boundaries[ZPLUS][1].type)
          {
            case DIRICHLET :
            case CONVECTIVE: bc1y[P-1][j][i] += coeffPlus*qy[P][j][i]; break;
            default        : break;
          }
        }
      }
    }
  }
  ierr = DMDAVecRestoreArray(vda, bc1yGlobal, &bc1y); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(vda, qyLocal, &qy); CHKERRQ(ierr);

  // W-FLUXES
  PetscReal ***bc1z;
  ierr = DMDAVecGetArray(wda, bc1zGlobal, &bc1z); CHKERRQ(ierr);
  PetscReal ***qz;
  ierr = DMDAVecGetArray(wda, qzLocal, &qz); CHKERRQ(ierr);
  ierr = DMDAGetCorners(wda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  ierr = DMDAGetInfo(wda, NULL, &M, &N, &P, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
  // x-faces
  if (flow->boundaries[XPLUS][2].type != PERIODIC) // don't update if the BC type is periodic
  {
    coeffMinus = alphaImplicit*nu*2.0/dxW[0]/(dxW[0]+dxW[1]);
    coeffPlus  = alphaImplicit*nu*2.0/dxW[M]/(dxW[M]+dxW[M-1]);
    // loop over all points on the x-face
    for (k=pstart; k<pstart+p; k++)
    {
      for (j=nstart; j<nstart+n; j++)
      {
        // -X
        if (mstart == 0) // if the -X face is in the current process
        {
          switch (flow->boundaries[XMINUS][2].type)
          {
            case DIRICHLET :
            case CONVECTIVE: bc1z[k][j][0] += coeffMinus*qz[k][j][-1]; break;
            default        : break;
          }
        }
        // +X
        if (mstart+m-1 == M-1) // if the +X face is in the current process
        {
          switch (flow->boundaries[XPLUS][2].type)
          {
            case DIRICHLET :
            case CONVECTIVE: bc1z[k][j][M-1] += coeffPlus*qz[k][j][M]; break;
            default        : break;
          }
        }
      }
    }
  }
  // y-faces
  if (flow->boundaries[YPLUS][2].type != PERIODIC) // don't update if the BC type is periodic
  {
    coeffMinus = alphaImplicit*nu*2.0/dyW[0]/(dyW[0]+dyW[1]);
    coeffPlus  = alphaImplicit*nu*2.0/dyW[N]/(dyW[N]+dyW[N-1]);
    // loop over all points on the y-face
    for (k=pstart; k<pstart+p; k++)
    {
      for (i=mstart; i<mstart+m; i++)
      {
        // -Y
        if (nstart == 0) // if the -Y face is in the current process
        {
          switch (flow->boundaries[YMINUS][2].type)
          {
            case DIRICHLET :
            case CONVECTIVE: bc1z[k][0][i] += coeffMinus*qz[k][-1][i]; break;
            default        : break;
          }
        }
        // +Y
        if (nstart+n-1 == N-1) // if the +Y face is in the current process
        {
          switch (flow->boundaries[YPLUS][2].type)
          {
            case DIRICHLET :
            case CONVECTIVE: bc1z[k][N-1][i] += coeffPlus*qz[k][N][i]; break;
            default        : break;
          }
        }
      }
    }
  }
  // z-faces
  if (flow->boundaries[ZPLUS][2].type != PERIODIC) // don't update if the BC type is periodic
  {
    coeffMinus = alphaImplicit*nu*2.0/dzW[0]/(dzW[0]+dzW[1]);
    coeffPlus  = alphaImplicit*nu*2.0/dzW[P]/(dzW[P]+dzW[P-1]);
    // loop over all points on the z-face
    for (j=nstart; j<nstart+n; j++)
    {
      for (i=mstart; i<mstart+m; i++)
      {
        // -Z
        if (pstart == 0) // if the -Z face is in the current process
        {
          switch (flow->boundaries[ZMINUS][2].type)
          {
            case DIRICHLET :
            case CONVECTIVE: bc1z[0][j][i] += coeffMinus*qz[-1][j][i]/(mesh->dx[i]*mesh->dy[j]); break;
            default        : break;
          }
        }
        // +Z
        if (pstart+p-1 == P-1) // if the +Z face is in the current process
        {
          switch (flow->boundaries[ZPLUS][2].type)
          {
            case DIRICHLET :
            case CONVECTIVE: bc1z[P-1][j][i] += coeffPlus*qz[P][j][i]/(mesh->dx[i]*mesh->dy[j]); break;
            default        : break;
          }
        }
      }
    }
  }
  ierr = DMDAVecRestoreArray(wda, bc1zGlobal, &bc1z); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(wda, qzLocal, &qz); CHKERRQ(ierr);

  ierr = DMCompositeRestoreAccess(qPack, bc1, &bc1xGlobal, &bc1yGlobal, &bc1zGlobal); CHKERRQ(ierr);

  return 0;
} // generateBC1