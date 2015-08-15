/***************************************************************************//**
 * \file updateBoundaryGhosts.inl
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the method `updateBoundaryGhosts` 
 *        of the class `NavierStokesSolver`.
 */


/**
 * \brief This function updates the values of the velocity on the domain boundaries.
 *
 * These values are stored in the ghost cells of the local vectors of the flux, 
 * but refer to the values at the locations of the boundaries. The following are 
 * the update rules for different boundary conditions:
 *
 * - **Dirichlet:** The values in the ghost cells are set to be the values of 
 *   the velocity at the boundary. Note that the velocity values are stored and 
 *   not the fluxes.
 * - **Neumann:** The values in the ghost cells are set to be the values of the 
 *   velocity at the grid points nearest to the boundary. Note that the velocity 
 *   values are stored and not the fluxes.
 * - **Convective:** The velocities near the boundary are convected outside 
 *   using the advection equation \f$ \frac{\partial u}{\partial t} + 
 *   u_\infty\frac{\partial u}{\partial x}=0 \f$. Using the discretized form of
 *   this equation, we can calculate the value of the velocity at the boundary.
 *   Note that the velocity values are stored and not the fluxes.
 * - **Periodic:** The values at the ghost cells are the fluxes from the points 
 *   nearest to the opposite edge of the domain. Here, the values do not 
 *   coincide with the periodic boundary, but instead are at the location of the
 *   grid point in a wrapped domain. This is handled automatically
 *   by PETSc when `DMCompositeScatter` is called, and so nothing explicit is
 *   done in this function.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::updateBoundaryGhosts()
{
  return 0;
} // updateBoundaryGhosts


// two-dimensional specialization
template <>
PetscErrorCode NavierStokesSolver<2>::updateBoundaryGhosts()
{
  PetscErrorCode ierr;
  
  PetscInt i, j,           // loop indices
           M, N,           // global number of values in each direction
           m, n,           // local number of values in each direction
           mstart, nstart; // starting indices of values on process
  PetscReal dt = parameters->dt,
            startStep = parameters->startStep;
  PetscReal beta; // convective speed

  // fluxes in x-direction
  PetscReal **qx;
  ierr = DMDAVecGetArray(uda, qxLocal, &qx); CHKERRQ(ierr);
  ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
  ierr = DMDAGetInfo(uda, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
  // left and right boundaries
  if (flow->boundaries[XMINUS][0].type != PERIODIC) // do not update if x-periodic
  {
    for (j=nstart; j<nstart+n; j++) // loop over y-direction
    {
      if (mstart == 0) // if left boundary on current process
      {
        switch (flow->boundaries[XMINUS][0].type)
        {
          case DIRICHLET:
            qx[j][-1] = flow->boundaries[XMINUS][0].value*mesh->dy[j];
            break;
          case CONVECTIVE: 
            beta = flow->boundaries[XMINUS][0].value * dt / mesh->dx[0];
            qx[j][-1] = (1-beta)*qx[j][-1] + beta*qx[j][0];
            if (timeStep == startStep)
              qx[j][-1] = qx[j][0];
            break;
          case NEUMANN: 
            qx[j][-1] = qx[j][0];
            break;
          default:
            break;
        }
      }
      if (mstart+m-1 == M-1) // if right boundary on current process
      {     
        switch (flow->boundaries[XPLUS][0].type)
        {         
          case DIRICHLET:
            qx[j][M] = flow->boundaries[XPLUS][0].value*mesh->dy[j];
            break;
          case CONVECTIVE:
            beta = flow->boundaries[XPLUS][0].value * dt / mesh->dx[M];
            qx[j][M] = (1-beta)*qx[j][M] + beta*qx[j][M-1];
            if (timeStep == startStep)
              qx[j][M] = qx[j][M-1];
            break;
          case NEUMANN:
            qx[j][M] = qx[j][M-1];
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
    for (i=mstart; i<mstart+m; i++) // loop over x-direction
    { 
      if (nstart == 0) // bottom boundary on current process
      {
        switch (flow->boundaries[YMINUS][0].type)
        {
          case DIRICHLET:
            qx[-1][i] = flow->boundaries[YMINUS][0].value;
            break;
          case CONVECTIVE:
            beta = flow->boundaries[YMINUS][1].value * dt / (0.5*mesh->dy[0]);
            qx[-1][i] = (1.0-beta)*qx[-1][i] + beta*qx[0][i]/mesh->dy[0];
            if (timeStep == startStep)
              qx[-1][i] = qx[0][i]/mesh->dy[0];
            break;
          case NEUMANN:
            qx[-1][i] = qx[0][i]/mesh->dy[0];
            break;
          default:
            break;
        }
      }
      if (nstart+n-1 == N-1) // top boundary on current process
      {
        switch (flow->boundaries[YPLUS][0].type)
        {
          case DIRICHLET:
            qx[N][i] = flow->boundaries[YPLUS][0].value;
            break;
          case CONVECTIVE:
            beta = flow->boundaries[YPLUS][1].value * dt / (0.5*mesh->dy[N]);
            qx[N][i] = (1.0-beta)*qx[N][i] + beta*qx[N-1][i]/mesh->dy[N-1];
            if (timeStep == startStep)
              qx[N][i] = qx[N-1][i]/mesh->dy[N-1];
            break;
          case NEUMANN:
            qx[N][i] = qx[N-1][i]/mesh->dy[N-1];
            break;
          default:
            break;
        }
      }
    }
  }
  ierr = DMDAVecRestoreArray(uda, qxLocal, &qx); CHKERRQ(ierr);
  
  // fluxes in y-direction
  PetscReal **qy;
  ierr = DMDAVecGetArray(vda, qyLocal, &qy); CHKERRQ(ierr);
  ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
  ierr = DMDAGetInfo(vda, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
  // left and right boundaries
  if (flow->boundaries[XMINUS][1].type != PERIODIC) // do not update if x-periodic
  {
    for (j=nstart; j<nstart+n; j++) // loop over y-direction
    {
      if (mstart == 0) // left boundary on current process
      {
        switch (flow->boundaries[XMINUS][1].type)
        {
          case DIRICHLET:
            qy[j][-1] = flow->boundaries[XMINUS][1].value;
            break;
          case CONVECTIVE: 
            beta = flow->boundaries[XMINUS][0].value * dt / (0.5*mesh->dx[0]);
            qy[j][-1] = (1.0-beta)*qy[j][-1] + beta*qy[j][0]/mesh->dx[0];
            if (timeStep == startStep)
              qy[j][-1] = qy[j][0]/mesh->dx[0];
            break;
          case NEUMANN:
            qy[j][-1] = qy[j][0]/mesh->dx[0];
            break;
          default:
            break;
        }
      }
      if (mstart+m-1 == M-1) // right boundary on current process
      {
        switch (flow->boundaries[XPLUS][1].type)
        {
          case DIRICHLET:
            qy[j][M] = flow->boundaries[XPLUS][1].value;
            break;
          case CONVECTIVE:
            beta = flow->boundaries[XPLUS][0].value * dt / (0.5*mesh->dx[M]);
            qy[j][M] = (1.0-beta)*qy[j][M] + beta*qy[j][M-1]/mesh->dx[M-1];
            if (timeStep == startStep)
              qy[j][M] = qy[j][M-1]/mesh->dx[M-1];
            break;
          case NEUMANN:
            qy[j][M] = qy[j][M-1]/mesh->dx[M-1];
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
    for (i=mstart; i<mstart+m; i++) // loop over x-direction
    { 
      if (nstart == 0) // bottom boundary on current process
      {
        switch (flow->boundaries[YMINUS][1].type)
        {
          case DIRICHLET:
            qy[-1][i] = flow->boundaries[YMINUS][1].value*mesh->dx[i];
            break;
          case CONVECTIVE:
            beta = flow->boundaries[YMINUS][1].value * dt / mesh->dy[0];
            qy[-1][i] = (1.0-beta)*qy[-1][i] + beta*qy[0][i];
            if (timeStep == startStep)
              qy[-1][i] = qy[0][i];
            break;
          case NEUMANN:
            qy[-1][i] = qy[0][i];
            break;
          default:
            break;
        }
      }
      if (nstart+n-1 == N-1) // top boundary on current process
      {
        switch (flow->boundaries[YPLUS][1].type)
        {
          case DIRICHLET:
            qy[N][i] = flow->boundaries[YPLUS][1].value*mesh->dx[i];
            break;
          case CONVECTIVE:
            beta = flow->boundaries[YPLUS][1].value * dt / mesh->dy[N];
            qy[N][i] = (1.0-beta)*qy[N][i] + beta*qy[N-1][i];
            if (timeStep == startStep)
              qy[N][i] = qy[N-1][i];
            break;
          case NEUMANN:
            qy[N][i] = qy[N-1][i];
            break;
          default:
            break;
        }
      }
    }
  }
  ierr = DMDAVecRestoreArray(vda, qyLocal, &qy); CHKERRQ(ierr);

  return 0;
} // updateBoundaryGhosts


// three-dimensional specialization
template <>
PetscErrorCode NavierStokesSolver<3>::updateBoundaryGhosts()
{
  PetscErrorCode ierr;

  PetscInt i, j, k,                // loop indices
           M, N, P,                // global number of values in each direction
           m, n, p,                // local number of values in each direction
           mstart, nstart, pstart; // stating indices
  
  PetscReal dt = parameters->dt,               // time-increment
            startStep = parameters->startStep; // starting time-step
  PetscReal beta; // convective speed

  // fluxes in x-direction
  PetscReal ***qx;
  ierr = DMDAVecGetArray(uda, qxLocal, &qx); CHKERRQ(ierr);
  ierr = DMDAGetCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  ierr = DMDAGetInfo(uda, NULL, &M, &N, &P, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
  // left and right boundaries
  if (flow->boundaries[XMINUS][0].type != PERIODIC) // do not update if x-periodic
  {
    for (k=pstart; k<pstart+p; k++) // loop over z-direction
    {
      for (j=nstart; j<nstart+n; j++) // loop over y-direction
      {
        if (mstart == 0) // left boundary on current process
        {
          switch (flow->boundaries[XMINUS][0].type)
          {
            case DIRICHLET:
              qx[k][j][-1] = flow->boundaries[XMINUS][0].value*(mesh->dy[j]*mesh->dz[k]);
              break;
            case CONVECTIVE:
              beta = flow->boundaries[XMINUS][0].value * dt / mesh->dx[0];
              qx[k][j][-1] = (1-beta)*qx[k][j][-1] + beta*qx[k][j][0];
              if (timeStep == startStep)
                qx[k][j][-1] = qx[k][j][0];
              break;
            case NEUMANN:
              qx[k][j][-1] = qx[k][j][0];
              break;
            default:
              break;
          }
        }
        if (mstart+m-1 == M-1) // right boundary on current process
        {
          switch (flow->boundaries[XPLUS][0].type)
          {
            case DIRICHLET:
              qx[k][j][M] = flow->boundaries[XPLUS][0].value*(mesh->dy[j]*mesh->dz[k]);
              break;
            case CONVECTIVE:
              beta = flow->boundaries[XPLUS][0].value * dt / mesh->dx[M];
              qx[k][j][M] = (1-beta)*qx[k][j][M] + beta*qx[k][j][M-1];
              if (timeStep == startStep)
                qx[k][j][M] = qx[k][j][M-1];
              break;
            case NEUMANN:
              qx[k][j][M] = qx[k][j][M-1];
              break;
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
    for (k=pstart; k<pstart+p; k++) // loop over z-direction
    {
      for (i=mstart; i<mstart+m; i++) // loop over x-direction
      {
        if (nstart == 0) // bottom boundary on current process
        {
          switch (flow->boundaries[YMINUS][0].type)
          {
            case DIRICHLET:
              qx[k][-1][i] = flow->boundaries[YMINUS][0].value;
              break;
            case CONVECTIVE:
              beta = flow->boundaries[YMINUS][1].value * dt / (0.5*mesh->dy[0]);
              qx[k][-1][i] = (1.0-beta)*qx[k][-1][i] + beta*qx[k][0][i]/(mesh->dy[0]*mesh->dz[k]);
              if (timeStep == startStep)
                qx[k][-1][i] = qx[k][0][i]/(mesh->dy[0]*mesh->dz[k]);
              break;
            case NEUMANN:
              qx[k][-1][i] = qx[k][0][i]/(mesh->dy[0]*mesh->dz[k]);
              break;
            default:
              break;
          }
        }
        if (nstart+n-1 == N-1) // top boundary on current process
        {
          switch (flow->boundaries[YPLUS][0].type)
          {
            case DIRICHLET:
              qx[k][N][i] = flow->boundaries[YPLUS][0].value;
              break;
            case CONVECTIVE:
              beta = flow->boundaries[YPLUS][1].value * dt / (0.5*mesh->dy[N]);
              qx[k][N][i] = (1.0-beta)*qx[k][N][i] + beta*qx[k][N-1][i]/(mesh->dy[N-1]*mesh->dz[k]);
              if (timeStep == startStep)
                qx[k][N][i] = qx[k][N-1][i]/(mesh->dy[N-1]*mesh->dz[k]);
              break;
            case NEUMANN:
              qx[k][N][i] = qx[k][N-1][i]/(mesh->dy[N-1]*mesh->dz[k]);
              break;
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
    for (j=nstart; j<nstart+n; j++) // loop over y-direction
    {
      for (i=mstart; i<mstart+m; i++) // loop over x-direction
      {
        if (pstart == 0) // back boundary on current process
        {
          switch (flow->boundaries[ZMINUS][0].type)
          {
            case DIRICHLET:
              qx[-1][j][i] = flow->boundaries[ZMINUS][0].value;
              break;
            case CONVECTIVE:
              beta = flow->boundaries[ZMINUS][2].value * dt / (0.5*mesh->dz[0]);
              qx[-1][j][i] = (1.0-beta)*qx[-1][j][i] + beta*qx[0][j][i]/(mesh->dy[j]*mesh->dz[0]);
              if (timeStep == startStep)
                qx[-1][j][i] = qx[0][j][i]/(mesh->dy[j]*mesh->dz[0]);
              break;
            case NEUMANN:
              qx[-1][j][i] = qx[0][j][i]/(mesh->dy[j]*mesh->dz[0]);
              break;
            default:
              break;
          }
        }
        if (pstart+p-1 == P-1) // fron boundary on current process
        {
          switch (flow->boundaries[ZPLUS][0].type)
          {
            case DIRICHLET:
              qx[P][j][i] = flow->boundaries[ZPLUS][0].value;
              break;
            case CONVECTIVE:
              beta = flow->boundaries[ZPLUS][2].value * dt / (0.5*mesh->dz[P]);
              qx[P][j][i] = (1.0-beta)*qx[P][j][i] + beta*qx[P-1][j][i]/(mesh->dy[j]*mesh->dz[P-1]);
              if (timeStep == startStep)
                qx[P][j][i] = qx[P-1][j][i]/(mesh->dy[j]*mesh->dz[P-1]);
              break;
            case NEUMANN:
              qx[P][j][i] = qx[P-1][j][i]/(mesh->dy[j]*mesh->dz[P-1]);
              break;
            default:
              break;
          }
        }
      }
    }
  }
  ierr = DMDAVecRestoreArray(uda, qxLocal, &qx); CHKERRQ(ierr);
  
  // fluxes in y-direction
  PetscReal ***qy;
  ierr = DMDAVecGetArray(vda, qyLocal, &qy); CHKERRQ(ierr);
  ierr = DMDAGetCorners(vda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  ierr = DMDAGetInfo(vda, NULL, &M, &N, &P, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
  // left and right boundaries
  if (flow->boundaries[XMINUS][1].type != PERIODIC) // do not update if x-periodic
  {
    for (k=pstart; k<pstart+p; k++) // loop over z-direction
    {
      for (j=nstart; j<nstart+n; j++) // loop over y-direction
      {
        if (mstart == 0) // left boundary on current process
        {
          switch (flow->boundaries[XMINUS][1].type)
          {
            case DIRICHLET:
              qy[k][j][-1] = flow->boundaries[XMINUS][1].value;
              break;
            case CONVECTIVE:
              beta = flow->boundaries[XMINUS][0].value * dt / (0.5*mesh->dx[0]);
              qy[k][j][-1] = (1.0-beta)*qy[k][j][-1] + beta*qy[k][j][0]/(mesh->dx[0]*mesh->dz[k]);
              if (timeStep == startStep)
                qy[k][j][-1] = qy[k][j][0]/(mesh->dx[0]*mesh->dz[k]);
              break;
            case NEUMANN:
              qy[k][j][-1] = qy[k][j][0]/(mesh->dx[0]*mesh->dz[k]);
              break;
            default:
              break;
          }
        }
        if (mstart+m-1 == M-1) // right boundary on current process
        {
          switch (flow->boundaries[XPLUS][1].type)
          {
            case DIRICHLET:
              qy[k][j][M] = flow->boundaries[XPLUS][1].value;
              break;
            case CONVECTIVE:
              beta = flow->boundaries[XPLUS][0].value * dt / (0.5*mesh->dx[M]);
              qy[k][j][M] = (1.0-beta)*qy[k][j][M] + beta*qy[k][j][M-1]/(mesh->dx[M-1]*mesh->dz[k]);
              if (timeStep == startStep)
                qy[k][j][M] = qy[k][j][M-1]/(mesh->dx[M-1]*mesh->dz[k]);
              break;
            case NEUMANN:
              qy[k][j][M] = qy[k][j][M-1]/(mesh->dx[M-1]*mesh->dz[k]);
              break;
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
    for (k=pstart; k<pstart+p; k++) // loop over z-direction
    {
      for (i=mstart; i<mstart+m; i++) // loop over x-direction
      {
        if (nstart == 0) // bottom boundary on current process
        {
          switch (flow->boundaries[YMINUS][1].type)
          {
            case DIRICHLET:
              qy[k][-1][i] = flow->boundaries[YMINUS][1].value*(mesh->dz[k]*mesh->dx[i]);
              break;
            case CONVECTIVE:
              beta = flow->boundaries[YMINUS][1].value * dt / mesh->dy[0];
              qy[k][-1][i] = (1.0-beta)*qy[k][-1][i] + beta*qy[k][0][i];
              if (timeStep == startStep)
                qy[k][-1][i] = qy[k][0][i];
              break;
            case NEUMANN:
              qy[k][-1][i] = qy[k][0][i];
              break;
            default:
              break;
          }
        }
        if (nstart+n-1 == N-1) // top boundary on current process
        {
          switch (flow->boundaries[YPLUS][1].type)
          {
            case DIRICHLET:
              qy[k][N][i] = flow->boundaries[YPLUS][1].value*(mesh->dz[k]*mesh->dx[i]);
              break;
            case CONVECTIVE:
              beta = flow->boundaries[YPLUS][1].value * dt / mesh->dy[N];
              qy[k][N][i] = (1.0-beta)*qy[k][N][i] + beta*qy[k][N-1][i];
              if (timeStep == startStep)
                qy[k][N][i] = qy[k][N-1][i];
              break;
            case NEUMANN:
              qy[k][N][i] = qy[k][N-1][i];
              break;
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
    for (j=nstart; j<nstart+n; j++) // loop over y-direction
    {
      for (i=mstart; i<mstart+m; i++) // loop over x-direction
      {
        if (pstart == 0) // back boundary on current process
        {
          switch (flow->boundaries[ZMINUS][1].type)
          {
            case DIRICHLET:
              qy[-1][j][i] = flow->boundaries[ZMINUS][1].value;
              break;
            case CONVECTIVE:
              beta = flow->boundaries[ZMINUS][2].value * dt / (0.5*mesh->dz[0]);
              qy[-1][j][i] = (1.0-beta)*qy[-1][j][i] + beta*qy[0][j][i]/(mesh->dx[i]*mesh->dz[0]);
              if (timeStep == startStep)
                qy[-1][j][i] = qy[0][j][i]/(mesh->dx[i]*mesh->dz[0]);
              break;
            case NEUMANN:
              qy[-1][j][i] = qy[0][j][i]/(mesh->dx[i]*mesh->dz[0]);
              break;
            default:
              break;
          }
        }
        if (pstart+p-1 == P-1) // front boundary on current process
        {
          switch (flow->boundaries[ZPLUS][1].type)
          {
            case DIRICHLET:
              qy[P][j][i] = flow->boundaries[ZPLUS][1].value;
              break;
            case CONVECTIVE:
              beta = flow->boundaries[ZPLUS][2].value * dt / (0.5*mesh->dz[P]);
              qy[P][j][i] = (1.0-beta)*qy[P][j][i] + beta*qy[P-1][j][i]/(mesh->dx[i]*mesh->dz[P-1]);
              if (timeStep == startStep)
                qy[P][j][i] = qy[P-1][j][i]/(mesh->dx[i]*mesh->dz[P-1]);
              break;
            case NEUMANN:
              qy[P][j][i] = qy[P-1][j][i]/(mesh->dx[i]*mesh->dz[P-1]);
              break;
            default:
              break;
          }
        }
      }
    }
  }
  ierr = DMDAVecRestoreArray(vda, qyLocal, &qy); CHKERRQ(ierr);

  // fluxes in z-direction
  PetscReal ***qz;
  ierr = DMDAVecGetArray(wda, qzLocal, &qz); CHKERRQ(ierr);
  ierr = DMDAGetCorners(wda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  ierr = DMDAGetInfo(wda, NULL, &M, &N, &P, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
  // left and right boundaries
  if (flow->boundaries[XMINUS][2].type != PERIODIC) // do not update if x-periodic
  {
    for (k=pstart; k<pstart+p; k++) // loop over z-direction
    {
      for (j=nstart; j<nstart+n; j++) // loop over y-direction
      {
        if (mstart == 0) // left boundary on current process
        {
          switch (flow->boundaries[XMINUS][2].type)
          {
            case DIRICHLET:
              qz[k][j][-1] = flow->boundaries[XMINUS][2].value;
              break;
            case CONVECTIVE:
              beta = flow->boundaries[XMINUS][0].value * dt / (0.5*mesh->dx[0]);
              qz[k][j][-1] = (1.0-beta)*qz[k][j][-1] + beta*qz[k][j][0]/(mesh->dx[0]*mesh->dy[j]);
              if (timeStep == startStep)
                qz[k][j][-1] = qz[k][j][0]/(mesh->dx[0]*mesh->dy[j]);
              break;
            case NEUMANN:
              qz[k][j][-1] = qz[k][j][0]/(mesh->dx[0]*mesh->dy[j]);
              break;
            default:
              break;
          }
        }
        if (mstart+m-1 == M-1) // right boundary on current process
        {
          switch (flow->boundaries[XPLUS][2].type)
          {
            case DIRICHLET:
              qz[k][j][M] = flow->boundaries[XPLUS][2].value;
              break;
            case CONVECTIVE:
              beta = flow->boundaries[XPLUS][0].value * dt / (0.5*mesh->dx[M]);
              qz[k][j][M] = (1.0-beta)*qz[k][j][M] + beta*qz[k][j][M-1]/(mesh->dx[M-1]*mesh->dy[j]);
              if (timeStep == startStep)
                qz[k][j][M] = qz[k][j][M-1]/(mesh->dx[M-1]*mesh->dy[j]);
              break;
            case NEUMANN:
              qz[k][j][M] = qz[k][j][M-1]/(mesh->dx[M-1]*mesh->dy[j]);
              break;
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
    for (k=pstart; k<pstart+p; k++) // loop over z-direction
    {
      for (i=mstart; i<mstart+m; i++) // loop over x-direction
      {
        if (nstart == 0) // bottom boundary on current process
        {
          switch (flow->boundaries[YMINUS][2].type)
          {
            case DIRICHLET:
              qz[k][-1][i] = flow->boundaries[YMINUS][2].value;
              break;
            case CONVECTIVE:
              beta = flow->boundaries[YMINUS][1].value * dt / (0.5*mesh->dy[0]);
              qz[k][-1][i] = (1.0-beta)*qz[k][-1][i] + beta*qz[k][0][i]/(mesh->dx[i]*mesh->dy[0]);
              if (timeStep == startStep)
                qz[k][-1][i] = qz[k][0][i]/(mesh->dx[i]*mesh->dy[0]);
              break;
            case NEUMANN:
              qz[k][-1][i] = qz[k][0][i]/(mesh->dx[i]*mesh->dy[0]);
              break;
            default:
              break;
          }
        }
        if (nstart+n-1 == N-1) // top boundary on current process
        {
          switch (flow->boundaries[YPLUS][2].type)
          {
            case DIRICHLET:
              qz[k][N][i] = flow->boundaries[YPLUS][2].value;
              break;
            case CONVECTIVE:
              beta = flow->boundaries[YPLUS][1].value * dt / (0.5*mesh->dy[N]);
              qz[k][N][i] = (1.0-beta)*qz[k][N][i] + beta*qz[k][N-1][i]/(mesh->dx[i]*mesh->dy[N-1]);
              if (timeStep == startStep)
                qz[k][N][i] = qz[k][N-1][i]/(mesh->dx[i]*mesh->dy[N-1]);
              break;
            case NEUMANN:
              qz[k][N][i] = qz[k][N-1][i]/(mesh->dx[i]*mesh->dy[N-1]);
              break;
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
    for (j=nstart; j<nstart+n; j++) // loop over y-direction
    {
      for (i=mstart; i<mstart+m; i++) // loop over x-direction
      {
        if (pstart == 0) // back boundary on current process
        {
          switch (flow->boundaries[ZMINUS][2].type)
          {
            case DIRICHLET:
              qz[-1][j][i] = flow->boundaries[ZMINUS][2].value*(mesh->dx[i]*mesh->dy[j]);
              break;
            case CONVECTIVE:
              beta = flow->boundaries[ZMINUS][2].value * dt / mesh->dz[0];
              qz[-1][j][i] = (1.0-beta)*qz[-1][j][i] + beta*qz[0][j][i];
              if (timeStep == startStep)
                qz[-1][j][i] = qz[0][j][i];
              break;
            case NEUMANN:
              qz[-1][j][i] = qz[0][j][i];
              break;
            default:
              break;
          }
        }
        if (pstart+p-1 == P-1) // front boundary on current process
        {
          switch (flow->boundaries[ZPLUS][2].type)
          {
            case DIRICHLET:
              qz[P][j][i] = flow->boundaries[ZPLUS][2].value*(mesh->dx[i]*mesh->dy[j]);
              break;
            case CONVECTIVE:
              beta = flow->boundaries[ZPLUS][2].value * dt / mesh->dz[P];
              qz[P][j][i] = (1.0-beta)*qz[P][j][i] + beta*qz[P-1][j][i];
              if (timeStep == startStep)
                qz[P][j][i] = qz[P-1][j][i];
              break;
            case NEUMANN:
              qz[P][j][i] = qz[P-1][j][i];
              break;
            default:
              break;
          }
        }
      }
    }
  }
  ierr = DMDAVecRestoreArray(wda, qzLocal, &qz); CHKERRQ(ierr);

  return 0;
} // updateBoundaryGhosts