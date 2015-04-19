/***************************************************************************//**
 * \file DiffusionTerm.cpp
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \brief Implementation of the methods of the class \c DiffusionTerm.
 */


#include "DiffusionTerm.h"

#include <sstream>

#include <petscdmcomposite.h>


/**
 * \brief Constructor - Initializes pointers.
 */
template <PetscInt dim>
DiffusionTerm<dim>::DiffusionTerm(std::string folder, 
                                  FlowDescription *FD, 
                                  SimulationParameters *SP, 
                                  CartesianMesh *CM) : NavierStokesSolver<dim>::NavierStokesSolver(folder, FD, SP, CM)
{
  rnExact = PETSC_NULL;
}

/**
 * \brief Initializes the solver for the diffusion equation.
 *
 * The explicit time-intregration coefficients for the convective terms 
 * are set to zero.
 */
template <PetscInt dim>
PetscErrorCode DiffusionTerm<dim>::initialize()
{
  PetscErrorCode ierr;

  ierr = NavierStokesSolver<dim>::initialize(); CHKERRQ(ierr);
  NavierStokesSolver<dim>::simParams->gamma = 0.0;
  NavierStokesSolver<dim>::simParams->zeta = 0.0;

  return ierr;
}

/**
 * \brief Initializes the fluxes with a sinusoidal solution.
 *
 * In terms of velocity, the initial conditions are:
 * \f[ u = \sin(\pi x)*\sin(\pi y)*\sin(\pi z) + 1.0 \f]
 * \f[ v = \sin(\pi x)*\sin(\pi y)*\sin(\pi z) \f]
 * \f[ w = \sin(\pi x)*\sin(\pi y)*\sin(\pi z) \f]
 * where (x, y, z) are the coordinates of the node.
 * To get the fluxes, the discrete velocities are multiplied by the area of 
 * the appropriate cell face.
 */
template <PetscInt dim>
PetscErrorCode DiffusionTerm<dim>::initializeFluxes()
{
  return 0;
}

template <>
PetscErrorCode DiffusionTerm<2>::initializeFluxes()
{
  PetscErrorCode ierr;

  PetscInt mBegin, nBegin,  // global indices of local lower-left corner
           m, n,            // number of local elements in each direction
           i, j;            // iteration indices
  PetscReal x, y;           // coordinates of a node

  Vec qxGlobal, qyGlobal;
  ierr = DMCompositeGetAccess(qPack, q, &qxGlobal, &qyGlobal); CHKERRQ(ierr);

  // fluxes in x-direction
  PetscReal **qx;
  ierr = DMDAVecGetArray(uda, qxGlobal, &qx); CHKERRQ(ierr);
  ierr = DMDAGetCorners(uda, &mBegin, &nBegin, NULL, &m, &n, NULL); CHKERRQ(ierr);
  PetscReal u;  // value of u-velocity
  for (j=nBegin; j<nBegin+n; j++)
  {
    for (i=mBegin; i<mBegin+m; i++)
    {
      x = mesh->x[i+1];
      y = 0.5*(mesh->y[j]+mesh->y[j+1]);
      u = sin(PETSC_PI*x)*sin(PETSC_PI*y) + 1.0;
      qx[j][i] = u*mesh->dy[j];
    }
  }
  ierr = DMDAVecRestoreArray(uda, qxGlobal, &qx); CHKERRQ(ierr);

  // fluxes in y-direction
  PetscReal **qy;
  ierr = DMDAVecGetArray(vda, qyGlobal, &qy); CHKERRQ(ierr);
  ierr = DMDAGetCorners(vda, &mBegin, &nBegin, NULL, &m, &n, NULL); CHKERRQ(ierr);
  PetscReal v;  // value of v-velocity
  for (j=nBegin; j<nBegin+n; j++)
  {
    for (i=mBegin; i<mBegin+m; i++)
    {
      x = 0.5*(mesh->x[i]+mesh->x[i+1]);
      y = mesh->y[j+1];
      v = sin(PETSC_PI*x)*sin(PETSC_PI*y);
      qy[j][i] = v*mesh->dx[i];
    }
  }
  ierr = DMDAVecRestoreArray(vda, qyGlobal, &qy); CHKERRQ(ierr);

  ierr = DMCompositeRestoreAccess(qPack, q, &qxGlobal, &qyGlobal); CHKERRQ(ierr);

  return ierr;
}

template <>
PetscErrorCode DiffusionTerm<3>::initializeFluxes()
{
  PetscErrorCode ierr;

  PetscInt mBegin, nBegin, pBegin,  // global indices of local lower-left corner
           m, n, p,                 // number of local elements in each direction
           i, j, k;                 // iteration indices
  PetscReal x, y, z;                // coordinates of a node

  Vec qxGlobal, qyGlobal, qzGlobal;
  ierr = DMCompositeGetAccess(qPack, q, &qxGlobal, &qyGlobal, &qzGlobal); CHKERRQ(ierr);

  // fluxes in x-direction
  PetscReal ***qx;
  ierr = DMDAVecGetArray(uda, qxGlobal, &qx); CHKERRQ(ierr);
  ierr = DMDAGetCorners(uda, &mBegin, &nBegin, &pBegin, &m, &n, &p); CHKERRQ(ierr);
  PetscReal u;  // value of u-velocity
  for (k=pBegin; k<pBegin+p; k++)
  {
    for (j=nBegin; j<nBegin+n; j++)
    {
      for (i=mBegin; i<mBegin+m; i++)
      {
        x = mesh->x[i+1];
        y = 0.5*(mesh->y[j]+mesh->y[j+1]);
        z = 0.5*(mesh->z[k]+mesh->z[k+1]);
        u = sin(PETSC_PI*x)*sin(PETSC_PI*y)*sin(PETSC_PI*z) + 1.0;
        qx[k][j][i] = u*mesh->dy[j]*mesh->dz[k];
      }
    }
  }
  ierr = DMDAVecRestoreArray(uda, qxGlobal, &qx); CHKERRQ(ierr);

  // fluxes in y-direction
  PetscReal ***qy;
  ierr = DMDAVecGetArray(vda, qyGlobal, &qy); CHKERRQ(ierr);
  ierr = DMDAGetCorners(vda, &mBegin, &nBegin, &pBegin, &m, &n, &p); CHKERRQ(ierr);
  PetscReal v;  // value of v-velocity
  for (k=pBegin; k<pBegin+p; k++)
  {
    for (j=nBegin; j<nBegin+n; j++)
    {
      for (i=mBegin; i<mBegin+m; i++)
      {
        x = 0.5*(mesh->x[i]+mesh->x[i+1]);
        y = mesh->y[j+1];
        z = 0.5*(mesh->z[k]+mesh->z[k+1]);
        v = sin(PETSC_PI*x)*sin(PETSC_PI*y)*sin(PETSC_PI*z);
        qy[k][j][i] = v*mesh->dx[i]*mesh->dz[k];
      }
    }
  }
  ierr = DMDAVecRestoreArray(vda, qyGlobal, &qy); CHKERRQ(ierr);

  // fluxes in z-direction
  PetscReal ***qz;
  ierr = DMDAVecGetArray(wda, qzGlobal, &qz); CHKERRQ(ierr);
  ierr = DMDAGetCorners(wda, &mBegin, &nBegin, &pBegin, &m, &n, &p); CHKERRQ(ierr);
  PetscReal w;  // value of w-velocity
  for (k=pBegin; k<pBegin+p; k++)
  {
    for (j=nBegin; j<nBegin+n; j++)
    {
      for (i=mBegin; i<mBegin+m; i++)
      {
        x = 0.5*(mesh->x[i]+mesh->x[i+1]);
        y = 0.5*(mesh->y[j]+mesh->y[j+1]);
        z = mesh->z[k+1];
        w = sin(PETSC_PI*x)*sin(PETSC_PI*y)*sin(PETSC_PI*z);
        qz[k][j][i] = w*mesh->dx[i]*mesh->dy[j];
      }
    }
  }
  ierr = DMDAVecRestoreArray(wda, qzGlobal, &qz); CHKERRQ(ierr);

  ierr = DMCompositeRestoreAccess(qPack, q, &qxGlobal, &qyGlobal, &qzGlobal); CHKERRQ(ierr);

  return ierr;
}

/**
 * \brief Computes the exact solution explicit diffusive terms.
 */
template <PetscInt dim>
PetscErrorCode DiffusionTerm<dim>::calculateExactSolution()
{
  return 0;
}

template <>
PetscErrorCode DiffusionTerm<2>::calculateExactSolution()
{
  PetscErrorCode ierr;
  
  PetscInt mBegin, nBegin,  // global indices of local lower-left node
           m, n,            // number of local elements in each direction
           i, j;            // iteration indices
  
  PetscReal x, y,           // coordinates of a node
            velocity,       // velocity value at the node
            diffusion;      // value of the diffusive term at the node
  
  PetscReal nu = flowDesc->nu,                // viscosity
            dt = simParams->dt,               // time-increment
            alpha = simParams->alphaExplicit; // explicit time-scheme coefficient for diffusive terms

  ierr = VecDuplicate(q, &rnExact); CHKERRQ(ierr);
  Vec rnExactXGlobal, rnExactYGlobal;
  ierr = DMCompositeGetAccess(qPack, rnExact, &rnExactXGlobal, &rnExactYGlobal); CHKERRQ(ierr);

  // exact explicit terms for u-nodes
  PetscReal **rnExactX;
  ierr = DMDAVecGetArray(uda, rnExactXGlobal, &rnExactX);
  ierr = DMDAGetCorners(uda, &mBegin, &nBegin, NULL, &m, &n, NULL); CHKERRQ(ierr);
  for (j=nBegin; j<nBegin+n; j++)
  {
    for (i=mBegin; i<mBegin+m; i++)
    {
      x = mesh->x[i+1];
      y = 0.5*(mesh->y[j]+mesh->y[j+1]);
      velocity = sin(PETSC_PI*x)*sin(PETSC_PI*y) + 1.0;
      diffusion = -2.0*nu*PETSC_PI*PETSC_PI*sin(PETSC_PI*x)*sin(PETSC_PI*y);
      rnExactX[j][i] = velocity/dt + alpha*diffusion;
    }
  }
  ierr = DMDAVecRestoreArray(uda, rnExactXGlobal, &rnExactX); CHKERRQ(ierr);

  // exact explicit terms for v-nodes
  PetscReal **rnExactY;
  ierr = DMDAVecGetArray(vda, rnExactYGlobal, &rnExactY);
  ierr = DMDAGetCorners(vda, &mBegin, &nBegin, NULL, &m, &n, NULL); CHKERRQ(ierr);
  for (j=nBegin; j<nBegin+n; j++)
  {
    for (i=mBegin; i<mBegin+m; i++)
    {
      x = 0.5*(mesh->x[i]+mesh->x[i+1]);
      y = mesh->y[j+1];
      velocity = sin(PETSC_PI*x)*sin(PETSC_PI*y);
      diffusion = -2.0*nu*PETSC_PI*PETSC_PI*sin(PETSC_PI*x)*sin(PETSC_PI*y);
      rnExactY[j][i] = velocity/dt + alpha*diffusion;
    }
  }
  ierr = DMDAVecRestoreArray(vda, rnExactYGlobal, &rnExactY); CHKERRQ(ierr);

  ierr = DMCompositeRestoreAccess(qPack, rnExact, &rnExactXGlobal, &rnExactYGlobal); CHKERRQ(ierr);

  return ierr;
}

template <>
PetscErrorCode DiffusionTerm<3>::calculateExactSolution()
{
  PetscErrorCode ierr;
  
  PetscInt mBegin, nBegin, pBegin,  // global indices of local lower-left node
           m, n, p,                 // number of local elements in each direction
           i, j, k;                 // iteration indices
  
  PetscReal x, y, z,                // coordinates of a node
            velocity,               // velocity value at the node
            diffusion;              // value of the diffusive term at the node
  
  PetscReal nu = flowDesc->nu,                // viscosity
            dt = simParams->dt,               // time-increment
            alpha = simParams->alphaExplicit; // explicit time-scheme coefficient for diffusive terms

  ierr = VecDuplicate(q, &rnExact); CHKERRQ(ierr);
  Vec rnExactXGlobal, rnExactYGlobal, rnExactZGlobal;
  ierr = DMCompositeGetAccess(qPack, rnExact, &rnExactXGlobal, &rnExactYGlobal, &rnExactZGlobal); CHKERRQ(ierr);

  // exact explicit terms for u-nodes
  PetscReal ***rnExactX;
  ierr = DMDAVecGetArray(uda, rnExactXGlobal, &rnExactX);
  ierr = DMDAGetCorners(uda, &mBegin, &nBegin, &pBegin, &m, &n, &p); CHKERRQ(ierr);
  for (k=pBegin; k<pBegin+p; k++)
  {
    for (j=nBegin; j<nBegin+n; j++)
    {
      for (i=mBegin; i<mBegin+m; i++)
      {
        x = mesh->x[i+1];
        y = 0.5*(mesh->y[j]+mesh->y[j+1]);
        z = 0.5*(mesh->z[k]+mesh->z[k+1]);
        velocity = sin(PETSC_PI*x)*sin(PETSC_PI*y)*sin(PETSC_PI*z) + 1.0;
        diffusion = -3.0*nu*PETSC_PI*PETSC_PI*sin(PETSC_PI*x)*sin(PETSC_PI*y)*sin(PETSC_PI*z);
        rnExactX[k][j][i] = velocity/dt + alpha*diffusion;
      }
    }
  }
  ierr = DMDAVecRestoreArray(uda, rnExactXGlobal, &rnExactX); CHKERRQ(ierr);

  // exact explicit terms for v-nodes
  PetscReal ***rnExactY;
  ierr = DMDAVecGetArray(vda, rnExactYGlobal, &rnExactY);
  ierr = DMDAGetCorners(vda, &mBegin, &nBegin, &pBegin, &m, &n, &p); CHKERRQ(ierr);
  for (k=pBegin; k<pBegin+p; k++)
  {
    for (j=nBegin; j<nBegin+n; j++)
    {
      for (i=mBegin; i<mBegin+m; i++)
      {
        x = 0.5*(mesh->x[i]+mesh->x[i+1]);
        y = mesh->y[j+1];
        z = 0.5*(mesh->z[k]+mesh->z[k+1]);
        velocity = sin(PETSC_PI*x)*sin(PETSC_PI*y)*sin(PETSC_PI*z);
        diffusion = -3.0*nu*PETSC_PI*PETSC_PI*sin(PETSC_PI*x)*sin(PETSC_PI*y)*sin(PETSC_PI*z);
        rnExactY[k][j][i] = velocity/dt + alpha*diffusion;
      }
    }
  }
  ierr = DMDAVecRestoreArray(vda, rnExactYGlobal, &rnExactY); CHKERRQ(ierr);

  // exact explicit terms for w-nodes
  PetscReal ***rnExactZ;
  ierr = DMDAVecGetArray(wda, rnExactZGlobal, &rnExactZ);
  ierr = DMDAGetCorners(wda, &mBegin, &nBegin, &pBegin, &m, &n, &p); CHKERRQ(ierr);
  for (k=pBegin; k<pBegin+p; k++)
  {
    for (j=nBegin; j<nBegin+n; j++)
    {
      for (i=mBegin; i<mBegin+m; i++)
      {
        x = 0.5*(mesh->x[i]+mesh->x[i+1]);
        y = 0.5*(mesh->y[j]+mesh->y[j+1]);
        z = mesh->z[k+1];
        velocity = sin(PETSC_PI*x)*sin(PETSC_PI*y)*sin(PETSC_PI*z);
        diffusion = -3.0*nu*PETSC_PI*PETSC_PI*sin(PETSC_PI*x)*sin(PETSC_PI*y)*sin(PETSC_PI*z);
        rnExactZ[k][j][i] = velocity/dt + alpha*diffusion;
      }
    }
  }
  ierr = DMDAVecRestoreArray(wda, rnExactZGlobal, &rnExactZ); CHKERRQ(ierr);

  ierr = DMCompositeRestoreAccess(qPack, rnExact, &rnExactXGlobal, &rnExactYGlobal, &rnExactZGlobal); CHKERRQ(ierr);

  return ierr;
}

/**
 * \brief Computes the relative L2 norm of the difference between the numerical
 *        and the exact solution of the explicit terms.
 */
template<PetscInt dim>
PetscErrorCode DiffusionTerm<dim>::calculateRelativeError()
{
  PetscErrorCode ierr;

  ierr = VecAXPY(NavierStokesSolver<dim>::rn, -1.0, rnExact);
  PetscReal l2NormDifference;
  ierr = VecNorm(NavierStokesSolver<dim>::rn, NORM_2, &l2NormDifference); CHKERRQ(ierr);
  PetscReal l2NormExact;
  ierr = VecNorm(rnExact, NORM_2, &l2NormExact); CHKERRQ(ierr);
  relativeError = l2NormDifference/l2NormExact;

  return ierr;
}

/**
 * \brief Writes the number of cells and the relative error in a file.
 *
 * In the file already exists, the values are appended; otherwise the file 
 * is created.
 */
template <PetscInt dim>
PetscErrorCode DiffusionTerm<dim>::writeRelativeError()
{
  PetscErrorCode ierr;

  PetscInt rank;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  if (rank == 0)
  {
    std::stringstream out;
    out << "./relative_errors_" << dim << "d.dat";
    std::ofstream fileStream;
    fileStream.open(out.str().c_str(), std::ios::out | std::ios::app);
    PetscInt nCells = NavierStokesSolver<dim>::mesh->nx*NavierStokesSolver<dim>::mesh->ny;
    if (dim == 3)
      nCells *= NavierStokesSolver<dim>::mesh->nz;
    fileStream << nCells << "\t" << relativeError << std::endl;
    fileStream.close();
  }

  return ierr;
}

/**
 * \brief Frees memory to avoid memory leaks
 */
template <PetscInt dim>
PetscErrorCode DiffusionTerm<dim>::finalize()
{
  PetscErrorCode ierr;

  if (rnExact != PETSC_NULL)
  {
    ierr = VecDestroy(&rnExact); CHKERRQ(ierr);
  }
  ierr = NavierStokesSolver<dim>::finalize(); CHKERRQ(ierr);
  
  return ierr;
}

template class DiffusionTerm<2>;
template class DiffusionTerm<3>;