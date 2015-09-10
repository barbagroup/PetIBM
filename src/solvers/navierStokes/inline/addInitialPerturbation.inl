/***************************************************************************//**
 * \file addInitialPerturbation.inl
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \brief Implementation of the method `addInitialPerturbation` 
 *        of the class `NavierStokesSolver`.
 */


/**
 * \brief Adds a perturbation to the initial flux field.
 *
 * This perturbation takes the form of a Taylor-Green vortex, whose amplitude
 * and frequency are stored in the object `flow`.
 *
 * The amplitude and the frequency of the vortices can be set in the input file 
 * `flowDescription.yaml` as following:
 * `initialPertubation: [amplitude, frequency]`
 * where `amplitude` and `frequency` are real values.
 *
 * In 2-D, this perturbation takes the form
 * \f[ u' = - amplitude \cos(X) \sin(Y) \f]
 * \f[ v' = amplitude \sin(X) \cos(Y) \f]
 *
 * and in 3-D
 * \f[ u' = - amplitude \cos(X) \sin(Y) \sin(Z) \f]
 * \f[ v' = amplitude \sin(X) \cos(Y) \sin(Z) \f]
 * \f[ w' = 0 \f]
 * 
 * where
 * \f[ X = 2\pi frequency (x-x_{start})/(x_{end}-x_{start}) \f]
 * \f[ X = 2\pi frequency (y-y_{start})/(y_{end}-y_{start}) \f]
 * \f[ X = 2\pi frequency (z-z_{start})/(z_{end}-z_{start}) \f]
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::addInitialPerturbation()
{
  return 0;
} // addInitialPerturbation


// two-dimensional specialization
template <>
PetscErrorCode NavierStokesSolver<2>::addInitialPerturbation()
{
  PetscErrorCode ierr;

  PetscInt i, j,           // loop indices
           m, n,           // local number of nodes along each direction
           mstart, nstart; // starting indices

  // Taylor-Green vortex perturbation
  PetscReal amplitude = flow->perturbationAmplitude,
            frequency = flow->perturbationFrequency;
  PetscReal X, Y; // scaled location

  // get access to individual packed vectors in their global representation
  Vec qxGlobal, qyGlobal;
  ierr = DMCompositeGetAccess(qPack, q, &qxGlobal, &qyGlobal); CHKERRQ(ierr);

  // fluxes in x-direction
  PetscReal **qx;
  ierr = DMDAVecGetArray(uda, qxGlobal, &qx); CHKERRQ(ierr);
  ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
  for (j=nstart; j<nstart+n; j++)
  {
    Y = 2.0*PETSC_PI*frequency * (0.5*(mesh->y[j]+mesh->y[j+1])-mesh->y[0]) / (mesh->y[mesh->ny]-mesh->y[0]);
    for (i=mstart; i<mstart+m; i++)
    {
      X = 2.0*PETSC_PI*frequency * (mesh->x[i+1]-mesh->x[0]) / (mesh->x[mesh->nx]-mesh->x[0]);    
      qx[j][i] += -amplitude*cos(X)*sin(Y) * mesh->dy[j];
    }
  }
  ierr = DMDAVecRestoreArray(uda, qxGlobal, &qx); CHKERRQ(ierr);
  // fluxes in y-direction
  PetscReal **qy;
  ierr = DMDAVecGetArray(vda, qyGlobal, &qy); CHKERRQ(ierr);
  ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
  for (j=nstart; j<nstart+n; j++)
  {
    Y = 2.0*PETSC_PI*frequency * (mesh->y[j]-mesh->y[0]) / (mesh->y[mesh->ny]-mesh->y[0]);
    for (i=mstart; i<mstart+m; i++)
    {
      X = 2.0*PETSC_PI*frequency * (0.5*(mesh->x[i]+mesh->x[i+1])-mesh->x[0]) / (mesh->x[mesh->nx]-mesh->x[0]);    
      qy[j][i] += amplitude*sin(X)*cos(Y) * mesh->dx[i];
    }
  }
  ierr = DMDAVecRestoreArray(vda, qyGlobal, &qy); CHKERRQ(ierr);

  ierr = DMCompositeRestoreAccess(qPack, q, &qxGlobal, &qyGlobal); CHKERRQ(ierr);

  return 0;
} // addInitialPerturbation


// three-dimensional specialization
template <>
PetscErrorCode NavierStokesSolver<3>::addInitialPerturbation()
{
  PetscErrorCode ierr;

  PetscInt i, j, k,                // loop indices
           m, n, p,                // local number of nodes along each direction
           mstart, nstart, pstart; // starting indices

  // Taylor-Green vortex perturbation
  PetscReal amplitude = flow->perturbationAmplitude,
            frequency = flow->perturbationFrequency;
  PetscReal X, Y, Z; // scaled location

  // get access to individual packed vectors in their global representation
  Vec qxGlobal, qyGlobal, qzGlobal;
  ierr = DMCompositeGetAccess(qPack, q, &qxGlobal, &qyGlobal, &qzGlobal); CHKERRQ(ierr);

  // fluxes in x-direction
  PetscReal ***qx;
  ierr = DMDAVecGetArray(uda, qxGlobal, &qx); CHKERRQ(ierr);
  ierr = DMDAGetCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  for (k=pstart; k<pstart+p; k++)
  {
    Z = 2.0*PETSC_PI*frequency * (0.5*(mesh->z[k]+mesh->z[k+1])-mesh->z[0]) / (mesh->z[mesh->nz]-mesh->z[0]);
    for (j=nstart; j<nstart+n; j++)
    {
      Y = 2.0*PETSC_PI*frequency * (0.5*(mesh->y[j]+mesh->y[j+1])-mesh->y[0]) / (mesh->y[mesh->ny]-mesh->y[0]);
      for (i=mstart; i<mstart+m; i++)
      {
        X = 2.0*PETSC_PI*frequency * (mesh->x[i+1]-mesh->x[0]) / (mesh->x[mesh->nx]-mesh->x[0]);    
        qx[k][j][i] += -amplitude*cos(X)*sin(Y)*sin(Z) * (mesh->dy[j]*mesh->dz[k]);
      }
    }
  }
  ierr = DMDAVecRestoreArray(uda, qxGlobal, &qx); CHKERRQ(ierr);
  // fluxes in y-direction
  PetscReal ***qy;
  ierr = DMDAVecGetArray(vda, qyGlobal, &qy); CHKERRQ(ierr);
  ierr = DMDAGetCorners(vda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  for (k=pstart; k<pstart+p; k++)
  {
    Z = 2.0*PETSC_PI*frequency * (0.5*(mesh->z[k]+mesh->z[k+1])-mesh->z[0]) / (mesh->z[mesh->nz]-mesh->z[0]);
    for (j=nstart; j<nstart+n; j++)
    {
      Y = 2.0*PETSC_PI*frequency * (mesh->y[j]-mesh->y[0]) / (mesh->y[mesh->ny]-mesh->y[0]);
      for (i=mstart; i<mstart+m; i++)
      {
        X = 2.0*PETSC_PI*frequency * (0.5*(mesh->x[i]+mesh->x[i+1])-mesh->x[0]) / (mesh->x[mesh->nx]-mesh->x[0]);    
        qy[k][j][i] += amplitude*sin(X)*cos(Y)*sin(Z) * (mesh->dx[i]*mesh->dz[k]);
      }
    }
  }
  ierr = DMDAVecRestoreArray(vda, qyGlobal, &qy); CHKERRQ(ierr);

  ierr = DMCompositeRestoreAccess(qPack, q, &qxGlobal, &qyGlobal, &qzGlobal); CHKERRQ(ierr);

  return 0;
} // addInitialPerturbation