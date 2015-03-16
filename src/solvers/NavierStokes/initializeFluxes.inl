/**
* Initialize the flux vector with the values of the velocity fluxes at the first
* time step. The initial values of the velocity from which the fluxes are 
* calculated are read from the FlowDescription object `flowDesc`. If specified,
* an initial perturbation is added to this velocity. In 2-D, this perturbation 
* takes the form
* \f[ u' = u_{perturbation} \cos\left(\frac{x}{2}\right) \sin(y) \f]
* \f[ v' = v_{perturbation} \sin(x) \cos\left(\frac{y}{2}\right) \f]
*
* and in 3-D:
* \f[ u' = u_{perturbation} \cos\left(\frac{x}{2}\right) \sin(y) \sin(z) \f]
* \f[ v' = v_{perturbation} \sin(x) \cos\left(\frac{y}{2}\right) \sin(z) \f]
* \f[ w' = w_{perturbation} \sin(x) \sin(y) \cos\left(\frac{z}{2}\right) \f]
*
* The above equations are calculated in the cube 
* \f$ [-\pi, \pi]\times[-\pi, \pi]\times[-\pi, \pi] \f$, and then mapped to the 
* cuboidal domain of the flow using a linear transformation:
* \f[ \left[-\pi,\pi\right] \rightarrow \left[x_{start}, x_{end}\right] \f]
* \f[ x \rightarrow -\pi + 2\pi\frac{x-x_{start}}{x_{end}-x_{start}} \f]
*
* Hence, the initial velocity of the flow field is given by:
* \f[ \vec{V}=\vec{V}_{initial}+\vec{V'} \f]
*
* When a simulation is restarted, the initial conditions are read from 
* previously saved data. The time step at which the data is read is specified 
* in the input file `simulationParameters.yaml`. The option `restart` is set to 
* `true`, and the option `startStep` specifies the time step from which the 
* simulation needs to be restarted.
*/
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::initializeFluxes()
{
  return 0;
}

template <>
PetscErrorCode NavierStokesSolver<2>::initializeFluxes()
{
  PetscErrorCode ierr;
  
  if(simParams->restart)
  {
    ierr = readFluxes(); CHKERRQ(ierr);
  }
  else
  {
    PetscInt  mstart, nstart, m, n;
    PetscReal **qx, **qy;
    Vec       qxGlobal, qyGlobal;
    PetscReal initVel[2]  = {flowDesc->initialVelocity[0], flowDesc->initialVelocity[1]};
    PetscReal initPert[2] = {flowDesc->initialPerturbation[0], flowDesc->initialPerturbation[1]};
    PetscReal width[2]    = {mesh->x[mesh->nx] - mesh->x[0], mesh->y[mesh->ny] - mesh->y[0]};

    ierr = DMCompositeGetAccess(qPack, q, &qxGlobal, &qyGlobal); CHKERRQ(ierr);

    // U-FLUXES
    ierr = DMDAVecGetArray(uda, qxGlobal, &qx); CHKERRQ(ierr);
    ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
    // Set interior values for u-fluxes
    for(PetscInt j=nstart; j<nstart+n; j++)
    {
      for(PetscInt i=mstart; i<mstart+m; i++)
      {
        PetscReal x = 2*PETSC_PI*(mesh->x[i+1] - mesh->x[0])/width[0],
                  y = 2*PETSC_PI*(0.5*(mesh->y[j]+mesh->y[j+1]) - mesh->y[0])/width[1];
        
        qx[j][i] = (initVel[0] + initPert[0]*sin(x)*cos(y))*mesh->dy[j];
      }
    }
    ierr = DMDAVecRestoreArray(uda, qxGlobal, &qx); CHKERRQ(ierr);
    
    // V-FLUXES
    ierr = DMDAVecGetArray(vda, qyGlobal, &qy); CHKERRQ(ierr);
    ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
    // Set interior values for v-fluxes
    for(PetscInt j=nstart; j<nstart+n; j++)
    {
      for(PetscInt i=mstart; i<mstart+m; i++)
      {
        PetscReal x = 2*PETSC_PI*(0.5*(mesh->x[i]+mesh->x[i+1]) - mesh->x[0])/width[0],
                  y = 2*PETSC_PI*(mesh->y[j+1] - mesh->y[0])/width[1];
        
        qy[j][i] = (initVel[1] + initPert[1]*cos(x)*sin(y))*mesh->dx[i];
      }
    }
    ierr = DMDAVecRestoreArray(vda, qyGlobal, &qy); CHKERRQ(ierr);

    ierr = DMCompositeRestoreAccess(qPack, q, &qxGlobal, &qyGlobal); CHKERRQ(ierr);
  }
  
  ierr = DMCompositeScatter(qPack, q, qxLocal, qyLocal); CHKERRQ(ierr);

  return 0;
}

template <>
PetscErrorCode NavierStokesSolver<3>::initializeFluxes()
{
  PetscErrorCode ierr;

  if(simParams->restart)
  {
    ierr = readFluxes(); CHKERRQ(ierr);
  }
  else
  {
    PetscInt  mstart, nstart, pstart, m, n, p;
    PetscReal ***qx, ***qy, ***qz;
    Vec       qxGlobal, qyGlobal, qzGlobal;
    PetscReal initVel[3]  = {flowDesc->initialVelocity[0], flowDesc->initialVelocity[1], flowDesc->initialVelocity[2]};
    PetscReal initPert[3] = {flowDesc->initialPerturbation[0], flowDesc->initialPerturbation[1], flowDesc->initialPerturbation[2]};
    PetscReal width[3]    = {mesh->x[mesh->nx] - mesh->x[0], mesh->y[mesh->ny] - mesh->y[0], mesh->z[mesh->nz] - mesh->z[0]};

    ierr = DMCompositeGetAccess(qPack, q, &qxGlobal, &qyGlobal, &qzGlobal); CHKERRQ(ierr);
    
    // U-FLUXES
    ierr = DMDAVecGetArray(uda, qxGlobal, &qx); CHKERRQ(ierr);
    ierr = DMDAGetCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
    // Set interior values for u-fluxes
    for(PetscInt k=pstart; k<pstart+p; k++)
    {
      for(PetscInt j=nstart; j<nstart+n; j++)
      {
        for(PetscInt i=mstart; i<mstart+m; i++)
        {
          PetscReal x = 2*PETSC_PI*(mesh->x[i+1] - mesh->x[0])/width[0],
                    y = 2*PETSC_PI*(0.5*(mesh->y[j]+mesh->y[j+1]) - mesh->y[0])/width[1],
                    z = 2*PETSC_PI*(0.5*(mesh->z[k]+mesh->z[k+1]) - mesh->z[0])/width[2];
          
          qx[k][j][i] = (initVel[0] + initPert[0]*sin(x)*cos(y)*cos(z))*(mesh->dy[j]*mesh->dz[k]);
        }
      }
    }
    ierr = DMDAVecRestoreArray(uda, qxGlobal, &qx); CHKERRQ(ierr);
    
    // V-FLUXES
    ierr = DMDAVecGetArray(vda, qyGlobal, &qy); CHKERRQ(ierr);
    ierr = DMDAGetCorners(vda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
    // Set interior values for v-fluxes
    for(PetscInt k=pstart; k<pstart+p; k++)
    {
      for(PetscInt j=nstart; j<nstart+n; j++)
      {
        for(PetscInt i=mstart; i<mstart+m; i++)
        {
          PetscReal x = 2*PETSC_PI*(0.5*(mesh->x[i]+mesh->x[i+1]) - mesh->x[0])/width[0],
                    y = 2*PETSC_PI*(mesh->y[j+1] - mesh->y[0])/width[1],
                    z = 2*PETSC_PI*(0.5*(mesh->z[k]+mesh->z[k+1]) - mesh->z[0])/width[2];
        
          qy[k][j][i] = (initVel[1] + initPert[1]*cos(x)*sin(y)*cos(z))*(mesh->dx[i]*mesh->dz[k]);
        }
      }
    }
    ierr = DMDAVecRestoreArray(vda, qyGlobal, &qy); CHKERRQ(ierr);

    // W-FLUXES
    ierr = DMDAVecGetArray(wda, qzGlobal, &qz); CHKERRQ(ierr);
    ierr = DMDAGetCorners(wda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
    // Set interior values for w-fluxes
    for(PetscInt k=pstart; k<pstart+p; k++)
    {
      for(PetscInt j=nstart; j<nstart+n; j++)
      {
        for(PetscInt i=mstart; i<mstart+m; i++)
        { 
          qz[k][j][i] = (initVel[2] + initPert[2]*0.0)*(mesh->dx[i]*mesh->dy[j]);
        }
      }
    }
    ierr = DMDAVecRestoreArray(wda, qzGlobal, &qz); CHKERRQ(ierr);

    ierr = DMCompositeRestoreAccess(qPack, q, &qxGlobal, &qyGlobal, &qzGlobal); CHKERRQ(ierr);
  }
  
  ierr = DMCompositeScatter(qPack, q, qxLocal, qyLocal, qzLocal); CHKERRQ(ierr);

  return 0;
}
