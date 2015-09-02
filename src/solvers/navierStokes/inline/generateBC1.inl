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
  PetscReal coeffMinus, coeffPlus;

  PetscReal *dx = &mesh->dx[0],
            *dy = &mesh->dy[0];

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
  if (flow->boundaries[XMINUS][0].type != PERIODIC) // do not update if x-periodic
  {
    coeffMinus = alpha*nu * 2.0/dx[0]/(dx[0]+dx[1]);
    coeffPlus = alpha*nu * 2.0/dx[M]/(dx[M-1]+dx[M]);
    for (j=nstart; j<nstart+n; j++) // loop over y-direction
    {
      if (mstart == 0) // left boundary on current process
      {
        switch (flow->boundaries[XMINUS][0].type)
        {
          case CONVECTIVE:
          case DIRICHLET:
            bc1x[j][0] += coeffMinus*qx[j][-1]/dy[j];
            break;
          default:
            break;
        }
      }
      if (mstart+m-1 == M-1) // right on current process
      {     
        switch (flow->boundaries[XPLUS][0].type)
        {
          case CONVECTIVE:
          case DIRICHLET:
            bc1x[j][M-1] += coeffPlus*qx[j][M]/dy[j];
            break;
          default:
            break;
        }
      }
    }
  }
  // bottom and top boundaries
  if (flow->boundaries[YMINUS][0].type != PERIODIC) // do not update if y-periodic
  {
    coeffMinus = alpha*nu * 2.0/(0.5*dy[0])/(0.5*dy[0]+0.5*(dy[0]+dy[1]));
    coeffPlus = alpha*nu * 2.0/(0.5*dy[N-1])/(0.5*(dy[N-2]+dy[N-1])+0.5*dy[N-1]);
    for (i=mstart; i<mstart+m; i++) // loop over x-direction
    { 
      if (nstart == 0) // bottom boundary on current boundary
      {
        switch (flow->boundaries[YMINUS][0].type)
        {
          case CONVECTIVE:
          case DIRICHLET:
            bc1x[0][i] += coeffMinus*qx[-1][i];
            break;
          default:
            break;
        }
      }
      if (nstart+n-1 == N-1) // top boundary on current process
      {
        switch (flow->boundaries[YPLUS][0].type)
        {
          case CONVECTIVE:
          case DIRICHLET:
            bc1x[N-1][i] += coeffPlus*qx[N][i];
            break;
          default:
            break;
        }
      }
    }
  }
  ierr = DMDAVecRestoreArray(uda, bc1xGlobal, &bc1x); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(uda, qxLocal, &qx); CHKERRQ(ierr);
  
  // fluxes in y-direction
  PetscReal **bc1y;
  ierr = DMDAVecGetArray(vda, bc1yGlobal, &bc1y); CHKERRQ(ierr);
  PetscReal **qy;
  ierr = DMDAVecGetArray(vda, qyLocal, &qy); CHKERRQ(ierr);
  ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
  ierr = DMDAGetInfo(vda, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
  // left and right boundaries
  if (flow->boundaries[XMINUS][1].type != PERIODIC) // do not update if x-periodic
  {
    coeffMinus = alpha*nu * 2.0/(0.5*dx[0])/(0.5*dx[0]+0.5*(dx[0]+dx[1]));
    coeffPlus = alpha*nu * 2.0/(0.5*dx[M-1])/(0.5*(dx[M-2]+dx[M-1])+0.5*dx[M-1]);
    for (j=nstart; j<nstart+n; j++) //loop over y-direction
    {
      if (mstart == 0) // left boundary on current process
      {
        switch (flow->boundaries[XMINUS][1].type)
        {
          case CONVECTIVE:
          case DIRICHLET:
            bc1y[j][0] += coeffMinus*qy[j][-1];
            break;
          default:
            break;
        }
      }
      if (mstart+m-1 == M-1) // right boundary on current process
      {
        switch (flow->boundaries[XPLUS][1].type)
        {
          case CONVECTIVE:
          case DIRICHLET:
            bc1y[j][M-1] += coeffPlus*qy[j][M];
            break;
          default:
            break;
        }
      }
    }
  }
  // bottom and top boundaries
  if (flow->boundaries[YMINUS][1].type != PERIODIC) // do not update if y-periodic
  {
    coeffMinus = alpha*nu * 2.0/dy[0]/(dy[0]+dy[1]);
    coeffPlus = alpha*nu * 2.0/dy[N]/(dy[N-1]+dy[N]);
    for (i=mstart; i<mstart+m; i++) // loop over x-direction
    { 
      if (nstart == 0) // bottom boundary on current process
      {
        switch (flow->boundaries[YMINUS][1].type)
        {
          case CONVECTIVE:
          case DIRICHLET:
            bc1y[0][i] += coeffMinus*qy[-1][i]/dx[i];
            break;
          default:
            break;
        }
      }
      if (nstart+n-1 == N-1) // top boundary on current boundary
      {
        switch (flow->boundaries[YPLUS][1].type)
        {
          case CONVECTIVE:
          case DIRICHLET:
            bc1y[N-1][i] += coeffPlus*qy[N][i]/dx[i];
            break;
          default:
            break;
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

  PetscInt i, j, k,                // loop indices 
           M, N, P,                // global number of nodes along each direction
           m, n, p,                // local number of nodes along each direction
           mstart, nstart, pstart; // starting indices
  
  PetscReal *dx = &mesh->dx[0],
            *dy = &mesh->dy[0],
            *dz = &mesh->dz[0];

  PetscReal nu = flow->nu; // viscosity
  PetscReal alpha = parameters->diffusion.coefficients[0]; // implicit diffusion coefficient
  PetscReal coeffMinus, coeffPlus;

  ierr = VecSet(bc1, 0.0); CHKERRQ(ierr);
  Vec bc1xGlobal, bc1yGlobal, bc1zGlobal;
  ierr = DMCompositeGetAccess(qPack, bc1, &bc1xGlobal, &bc1yGlobal, &bc1zGlobal); CHKERRQ(ierr);
  
  // fluxes in x-direction
  PetscReal ***bc1x;
  ierr = DMDAVecGetArray(uda, bc1xGlobal, &bc1x); CHKERRQ(ierr);
  PetscReal ***qx;
  ierr = DMDAVecGetArray(uda, qxLocal, &qx); CHKERRQ(ierr);
  ierr = DMDAGetCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  ierr = DMDAGetInfo(uda, NULL, &M, &N, &P, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
  // left and right boundaries
  if (flow->boundaries[XMINUS][0].type != PERIODIC) // do not update if x-periodic
  {
    coeffMinus = alpha*nu * 2.0/dx[0]/(dx[0]+dx[1]);
    coeffPlus = alpha*nu * 2.0/dx[M]/(dx[M-1]+dx[M]);
    for (k=pstart; k<pstart+p; k++) // loop over z-direction
    {
      for (j=nstart; j<nstart+n; j++) // loop over y-direction
      {
        if (mstart == 0) // left boundary on current process
        {
          switch (flow->boundaries[XMINUS][0].type)
          {
            case DIRICHLET:
              bc1x[k][j][0] += coeffMinus*qx[k][j][-1]/(dy[j]*dz[k]);
              break;
            case CONVECTIVE:
            default:
              break;
          }
        }
        if (mstart+m-1 == M-1) // right boundary on current process
        {
          switch (flow->boundaries[XPLUS][0].type)
          {
            case DIRICHLET:
              bc1x[k][j][M-1] += coeffPlus*qx[k][j][M]/(dy[j]*dz[k]);
              break;
            case CONVECTIVE:
            default:
              break;
          }
        }
      }
    }
  }
  // bottom and top boundaries
  if (flow->boundaries[YMINUS][0].type != PERIODIC) // do not update if y-periodic
  {
    coeffMinus = alpha*nu * 2.0/(0.5*dy[0])/(0.5*dy[0]+0.5*(dy[0]+dy[1]));
    coeffPlus = alpha*nu * 2.0/(0.5*dy[N])/(0.5*(dy[N-1]+dy[N])+0.5*dy[N]);
    for (k=pstart; k<pstart+p; k++) // loop over z-direction
    {
      for (i=mstart; i<mstart+m; i++) // loop over y-direction
      {
        if (nstart == 0) // bottom boundary on current process
        {
          switch (flow->boundaries[YMINUS][0].type)
          {
            case DIRICHLET:
              bc1x[k][0][i] += coeffMinus*qx[k][-1][i];
              break;
            case CONVECTIVE:
            default:
              break;
          }
        }
        if (nstart+n-1 == N-1) // top boundary on current process
        {
          switch (flow->boundaries[YPLUS][0].type)
          {
            case DIRICHLET:
              bc1x[k][N-1][i] += coeffPlus*qx[k][N][i];
              break;
            case CONVECTIVE:
            default:
              break;
          }
        }
      }
    }
  }
  // back and front boundaries
  if (flow->boundaries[ZMINUS][0].type != PERIODIC) // do not update if z-periodic
  {
    coeffMinus = alpha*nu * 2.0/(0.5*dz[0])/(0.5*dz[0]+0.5*(dz[0]+dz[1]));
    coeffPlus = alpha*nu * 2.0/(0.5*dz[P])/(0.5*(dz[P-1]+dz[P])+0.5*dz[P]);
    for (j=nstart; j<nstart+n; j++) // loop over y-direction
    {
      for (i=mstart; i<mstart+m; i++) // loop over x-direction
      {
        if (pstart == 0) // back boundary on current process
        {
          switch (flow->boundaries[ZMINUS][0].type)
          {
            case DIRICHLET:
              bc1x[0][j][i] += coeffMinus*qx[-1][j][i];
              break;
            case CONVECTIVE:
            default:
              break;
          }
        }
        if (pstart+p-1 == P-1) // front boundary on current process
        {
          switch (flow->boundaries[ZPLUS][0].type)
          {
            case DIRICHLET:
              bc1x[P-1][j][i] += coeffPlus*qx[P][j][i];
              break;
            case CONVECTIVE:
            default:
              break;
          }
        }
      }
    }
  }
  ierr = DMDAVecRestoreArray(uda, bc1xGlobal, &bc1x); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(uda, qxLocal, &qx); CHKERRQ(ierr);
  
  // fluxes in y-direction
  PetscReal ***bc1y;
  ierr = DMDAVecGetArray(vda, bc1yGlobal, &bc1y); CHKERRQ(ierr);
  PetscReal ***qy;
  ierr = DMDAVecGetArray(vda, qyLocal, &qy); CHKERRQ(ierr);
  ierr = DMDAGetCorners(vda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  ierr = DMDAGetInfo(vda, NULL, &M, &N, &P, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
  // left and right boundaries
  if (flow->boundaries[XMINUS][1].type != PERIODIC) // do not update if x-periodic
  {
    coeffMinus = alpha*nu * 2.0/(0.5*dx[0])/(0.5*dx[0]+0.5*(dx[0]+dx[1]));
    coeffPlus = alpha*nu * 2.0/(0.5*dx[M])/(0.5*(dx[M-1]+dx[M])+0.5*dx[M]);
    for (k=pstart; k<pstart+p; k++) // loop over z-direction
    {
      for (j=nstart; j<nstart+n; j++) // loop over y-direction
      {
        if (mstart == 0) // right boundary on current process
        {
          switch (flow->boundaries[XMINUS][1].type)
          {
            case DIRICHLET:
              bc1y[k][j][0] += coeffMinus*qy[k][j][-1];
              break;
            case CONVECTIVE:
            default:
              break;
          }
        }
        if (mstart+m-1 == M-1) // right boundary on current process
        {
          switch (flow->boundaries[XPLUS][1].type)
          {
            case DIRICHLET:
              bc1y[k][j][M-1] += coeffPlus*qy[k][j][M];
              break;
            case CONVECTIVE:
            default:
              break;
          }
        }
      }
    }
  }
  // bottom and top boundaries
  if (flow->boundaries[YMINUS][1].type != PERIODIC) // do not update if y-periodic
  {
    coeffMinus = alpha*nu * 2.0/dy[0]/(dy[0]+dy[1]);
    coeffPlus = alpha*nu * 2.0/dy[N]/(dy[N-1]+dy[N]);
    for (k=pstart; k<pstart+p; k++) // loop over z-direction
    {
      for (i=mstart; i<mstart+m; i++) // loop over x-direction
      {
        if (nstart == 0) // bottom boundary on current process
        {
          switch (flow->boundaries[YMINUS][1].type)
          {
            case DIRICHLET:
              bc1y[k][0][i] += coeffMinus*qy[k][-1][i]/(dx[i]*dz[k]);
              break;
            case CONVECTIVE:
            default:
              break;
          }
        }
        if (nstart+n-1 == N-1) // top boundary on current process
        {
          switch (flow->boundaries[YPLUS][1].type)
          {
            case DIRICHLET:
              bc1y[k][N-1][i] += coeffPlus*qy[k][N][i]/(dx[i]*dz[k]);
              break;
            case CONVECTIVE:
            default:
              break;
          }
        }
      }
    }
  }
  // back and front boundaries
  if (flow->boundaries[ZMINUS][1].type != PERIODIC) // do not update if z-periodic
  {
    coeffMinus = alpha*nu * 2.0/(0.5*dz[0])/(0.5*dz[0]+0.5*(dz[0]+dz[1]));
    coeffPlus = alpha*nu * 2.0/(0.5*dz[P])/(0.5*(dz[P-1]+dz[P])+0.5*dz[P]);
    for (j=nstart; j<nstart+n; j++) // loop over y-direction
    {
      for (i=mstart; i<mstart+m; i++) // loop over x-direction
      {
        if (pstart == 0) // back boundary on current process
        {
          switch (flow->boundaries[ZMINUS][1].type)
          {
            case DIRICHLET:
              bc1y[0][j][i] += coeffMinus*qy[-1][j][i];
              break;
            case CONVECTIVE:
            default:
              break;
          }
        }
        if (pstart+p-1 == P-1) // front boundary on current process
        {
          switch (flow->boundaries[ZPLUS][1].type)
          {
            case DIRICHLET:
              bc1y[P-1][j][i] += coeffPlus*qy[P][j][i];
              break;
            case CONVECTIVE:
            default:
              break;
          }
        }
      }
    }
  }
  ierr = DMDAVecRestoreArray(vda, bc1yGlobal, &bc1y); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(vda, qyLocal, &qy); CHKERRQ(ierr);

  // fluxes in z-direction
  PetscReal ***bc1z;
  ierr = DMDAVecGetArray(wda, bc1zGlobal, &bc1z); CHKERRQ(ierr);
  PetscReal ***qz;
  ierr = DMDAVecGetArray(wda, qzLocal, &qz); CHKERRQ(ierr);
  ierr = DMDAGetCorners(wda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  ierr = DMDAGetInfo(wda, NULL, &M, &N, &P, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
  // left and right boundaries
  if (flow->boundaries[XMINUS][2].type != PERIODIC) // do not update if x-periodic
  {
    coeffMinus = alpha*nu * 2.0/(0.5*dx[0])/(0.5*dx[0]+0.5*(dx[0]+dx[1]));
    coeffPlus = alpha*nu * 2.0/(0.5*dx[M])/(0.5*(dx[M-1]+dx[M])+0.5*dx[M]);
    for (k=pstart; k<pstart+p; k++) // loop over z-direction
    {
      for (j=nstart; j<nstart+n; j++) // loop over y-direction
      {
        if (mstart == 0) // left boundary on current process
        {
          switch (flow->boundaries[XMINUS][2].type)
          {
            case DIRICHLET:
              bc1z[k][j][0] += coeffMinus*qz[k][j][-1];
              break;
            case CONVECTIVE:
            default:
              break;
          }
        }
        if (mstart+m-1 == M-1) // right boundary on current process
        {
          switch (flow->boundaries[XPLUS][2].type)
          {
            case DIRICHLET:
              bc1z[k][j][M-1] += coeffPlus*qz[k][j][M];
              break;
            case CONVECTIVE:
            default:
              break;
          }
        }
      }
    }
  }
  // bottom and top boundaries
  if (flow->boundaries[YMINUS][2].type != PERIODIC) // do not update if y-periodic
  {
    coeffMinus = alpha*nu * 2.0/(0.5*dy[0])/(0.5*dy[0]+0.5*(dy[0]+dy[1]));
    coeffPlus = alpha*nu * 2.0/(0.5*dy[N])/(0.5*(dy[N-1]+dy[N])+0.5*dy[N]);
    for (k=pstart; k<pstart+p; k++) // loop over z-direction
    {
      for (i=mstart; i<mstart+m; i++) // loop over x-direction
      {
        if (nstart == 0) // bottom boundary on current process
        {
          switch (flow->boundaries[YMINUS][2].type)
          {
            case DIRICHLET:
              bc1z[k][0][i] += coeffMinus*qz[k][-1][i];
              break;
            case CONVECTIVE:
            default:
              break;
          }
        }
        if (nstart+n-1 == N-1) // top boundary on current process
        {
          switch (flow->boundaries[YPLUS][2].type)
          {
            case DIRICHLET:
              bc1z[k][N-1][i] += coeffPlus*qz[k][N][i];
              break;
            case CONVECTIVE:
            default:
              break;
          }
        }
      }
    }
  }
  // back and front boundaries
  if (flow->boundaries[ZMINUS][2].type != PERIODIC) // do not update if z-periodic
  {
    coeffMinus = alpha*nu * 2.0/dz[0]/(dz[0]+dz[1]);
    coeffPlus = alpha*nu * 2.0/dz[P]/(dz[P-1]+dz[P]);
    for (j=nstart; j<nstart+n; j++) // loop over y-direction
    {
      for (i=mstart; i<mstart+m; i++) // loop over x-direction
      {
        if (pstart == 0) // back boundary on current process
        {
          switch (flow->boundaries[ZMINUS][2].type)
          {
            case DIRICHLET:
              bc1z[0][j][i] += coeffMinus*qz[-1][j][i]/(dx[i]*dy[j]);
              break;
            case CONVECTIVE:
            default:
              break;
          }
        }
        if (pstart+p-1 == P-1) // front boundary on current process
        {
          switch (flow->boundaries[ZPLUS][2].type)
          {
            case DIRICHLET:
              bc1z[P-1][j][i] += coeffPlus*qz[P][j][i]/(dx[i]*dy[j]);
              break;
            case CONVECTIVE:
            default:
              break;
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