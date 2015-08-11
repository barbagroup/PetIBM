/***************************************************************************//**
 * \file initializeFluxes.inl
 * \author Anush Kirshnan (anush@bu.edu), Olivier Mesnard (mesnardo@gwu.edu)
 * \brief Implementation of the method `initializeFluxes`.
 */


/**
 * \brief Computes the initial fluxes with an optional initial perturbation.
 * 
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
 * In the input file `simulationParameters.yaml`, if `startStep` is higher than 
 * `0`, the simulation is restarted, reading the numerical solution at the 
 * time-step specified.
 *
 * In the input file `flowDescription.yaml`, if the parameter `initialCustomField` 
 * is set to `true`, the initial field is read from the solution folder `0000000`.
 * However if the input `startStep` is higher than `0`, the simulation is considered
 * as restarted.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::initializeFluxes()
{
  return 0;
} // initializeFluxes


// specialization in 2 dimensions
template <>
PetscErrorCode NavierStokesSolver<2>::initializeFluxes()
{
  PetscErrorCode ierr;
  
  if (parameters->startStep > 0 || flow->initialCustomField)
  {
    ierr = readFluxes(); CHKERRQ(ierr);
  }
  else
  { 
    PetscReal initialVelocity[2]  = {flow->initialVelocity[0], flow->initialVelocity[1]};
    PetscReal width[2] = {mesh->x[mesh->nx] - mesh->x[0], mesh->y[mesh->ny] - mesh->y[0]};

    // Taylor-Green vortex perturbation
    PetscReal amplitude = flow->perturbationAmplitude,
              frequency = flow->perturbationFrequency;
    // size of the Taylor-Green vortex
    PetscReal X1 = 0.0,
              X2 = 2.0*PETSC_PI,
              Y1 = 0.0,
              Y2 = 2.0*PETSC_PI;

    PetscInt mstart, nstart, m, n;

    Vec qxGlobal, qyGlobal;
    ierr = DMCompositeGetAccess(qPack, q, &qxGlobal, &qyGlobal); CHKERRQ(ierr);

    // fluxes in x-direction
    PetscReal **qx;
    ierr = DMDAVecGetArray(uda, qxGlobal, &qx); CHKERRQ(ierr);
    ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
    // set interior values of fluxes in x-direction
    for (PetscInt j=nstart; j<nstart+n; j++)
    {
      for (PetscInt i=mstart; i<mstart+m; i++)
      {
        PetscReal x = X1 + (X2-X1)*frequency*(mesh->x[i+1] - mesh->x[0])/width[0],
                  y = Y1 + (Y2-Y1)*frequency*(0.5*(mesh->y[j]+mesh->y[j+1]) - mesh->y[0])/width[1];
        qx[j][i] = (initialVelocity[0] - amplitude*cos(x)*sin(y))*mesh->dy[j];
      }
    }
    ierr = DMDAVecRestoreArray(uda, qxGlobal, &qx); CHKERRQ(ierr);
    
    // fluxes in y-direction
    PetscReal **qy;
    ierr = DMDAVecGetArray(vda, qyGlobal, &qy); CHKERRQ(ierr);
    ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
    // set interior values of fluxes in y-direction
    for (PetscInt j=nstart; j<nstart+n; j++)
    {
      for (PetscInt i=mstart; i<mstart+m; i++)
      {
        PetscReal x = X1 + (X2-X1)*frequency*(0.5*(mesh->x[i]+mesh->x[i+1]) - mesh->x[0])/width[0],
                  y = Y1 + (Y2-Y1)*frequency*(mesh->y[j+1] - mesh->y[0])/width[1];
        qy[j][i] = (initialVelocity[1] + amplitude*sin(x)*cos(y))*mesh->dx[i];
      }
    }
    ierr = DMDAVecRestoreArray(vda, qyGlobal, &qy); CHKERRQ(ierr);

    ierr = DMCompositeRestoreAccess(qPack, q, &qxGlobal, &qyGlobal); CHKERRQ(ierr);
  }

  ierr = DMCompositeScatter(qPack, q, qxLocal, qyLocal); CHKERRQ(ierr);

  return 0;
} // initializeFluxes


// specialization in 3 dimensions
template <>
PetscErrorCode NavierStokesSolver<3>::initializeFluxes()
{
  PetscErrorCode ierr;

  if (parameters->startStep > 0 || flow->initialCustomField)
  {
    ierr = readFluxes(); CHKERRQ(ierr);
  }
  else
  {
    PetscReal initialVelocity[3]  = {flow->initialVelocity[0], flow->initialVelocity[1], flow->initialVelocity[2]};
    PetscReal width[3] = {mesh->x[mesh->nx] - mesh->x[0], mesh->y[mesh->ny] - mesh->y[0], mesh->z[mesh->nz] - mesh->z[0]};

    // Taylor-Green vortex perturbation
    PetscReal amplitude = flow->perturbationAmplitude,
              frequency = flow->perturbationFrequency;
    PetscReal X1 = 0.0,
              X2 = 2.0*PETSC_PI,
              Y1 = 0.0,
              Y2 = 2.0*PETSC_PI,
              Z1 = 0.0,
              Z2 = 2.0*PETSC_PI;

    PetscInt  mstart, nstart, pstart, m, n, p;

    Vec qxGlobal, qyGlobal, qzGlobal;
    ierr = DMCompositeGetAccess(qPack, q, &qxGlobal, &qyGlobal, &qzGlobal); CHKERRQ(ierr);
    
    // fluxes in x-direction
    PetscReal ***qx;
    ierr = DMDAVecGetArray(uda, qxGlobal, &qx); CHKERRQ(ierr);
    ierr = DMDAGetCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
    // set interior values of fluxes in x-direction
    for (PetscInt k=pstart; k<pstart+p; k++)
    {
      for (PetscInt j=nstart; j<nstart+n; j++)
      {
        for (PetscInt i=mstart; i<mstart+m; i++)
        {
          PetscReal x = X1 + (X2-X1)*frequency*(mesh->x[i+1] - mesh->x[0])/width[0],
                    y = Y1 + (Y2-Y1)*frequency*(0.5*(mesh->y[j]+mesh->y[j+1]) - mesh->y[0])/width[1],
                    z = Z1 + (Z2-Z1)*frequency*(0.5*(mesh->z[k]+mesh->z[k+1]) - mesh->z[0])/width[2];     
          qx[k][j][i] = (initialVelocity[0] - amplitude*cos(x)*sin(y)*sin(z))*(mesh->dy[j]*mesh->dz[k]);
        }
      }
    }
    ierr = DMDAVecRestoreArray(uda, qxGlobal, &qx); CHKERRQ(ierr);
    
    // fluxes in y-direction
    PetscReal ***qy;
    ierr = DMDAVecGetArray(vda, qyGlobal, &qy); CHKERRQ(ierr);
    ierr = DMDAGetCorners(vda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
    // set interior values of fluxes in y-direction
    for (PetscInt k=pstart; k<pstart+p; k++)
    {
      for (PetscInt j=nstart; j<nstart+n; j++)
      {
        for (PetscInt i=mstart; i<mstart+m; i++)
        {
          PetscReal x = X1 + (X2-X1)*frequency*(0.5*(mesh->x[i]+mesh->x[i+1]) - mesh->x[0])/width[0],
                    y = Y1 + (Y2-Y1)*frequency*(mesh->y[j+1] - mesh->y[0])/width[1],
                    z = Z1 + (Z2-Z1)*frequency*(0.5*(mesh->z[k]+mesh->z[k+1]) - mesh->z[0])/width[2];
          qy[k][j][i] = (initialVelocity[1] + amplitude*sin(x)*cos(y)*sin(z))*(mesh->dx[i]*mesh->dz[k]);
        }
      }
    }
    ierr = DMDAVecRestoreArray(vda, qyGlobal, &qy); CHKERRQ(ierr);

    // fluxes in z-direction
    PetscReal ***qz;
    ierr = DMDAVecGetArray(wda, qzGlobal, &qz); CHKERRQ(ierr);
    ierr = DMDAGetCorners(wda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
    // set interior values of fluxes in z-direction
    for (PetscInt k=pstart; k<pstart+p; k++)
    {
      for (PetscInt j=nstart; j<nstart+n; j++)
      {
        for (PetscInt i=mstart; i<mstart+m; i++)
        { 
          qz[k][j][i] = initialVelocity[2]*(mesh->dx[i]*mesh->dy[j]);
        }
      }
    }
    ierr = DMDAVecRestoreArray(wda, qzGlobal, &qz); CHKERRQ(ierr);

    ierr = DMCompositeRestoreAccess(qPack, q, &qxGlobal, &qyGlobal, &qzGlobal); CHKERRQ(ierr);
  }
  
  ierr = DMCompositeScatter(qPack, q, qxLocal, qyLocal, qzLocal); CHKERRQ(ierr);

  return 0;
} // initializeFluxes