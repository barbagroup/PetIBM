/***************************************************************************//**
 * \file ConvectiveTerm.cpp
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \brief Implementation of the methods of the class \c ConvectiveTerm.
 */


#include "ConvectiveTerm.h"

#include <sstream>

#include <petscdmcomposite.h>


/**
 * \brief Constructor - Initializes pointers.
 */
template <PetscInt dim>
ConvectiveTerm<dim>::ConvectiveTerm(CartesianMesh *cartesianMesh, 
                                    FlowDescription<dim> *flowDescription, 
                                    SimulationParameters *simulationParameters) 
                   : NavierStokesSolver<dim>::NavierStokesSolver(cartesianMesh, 
                                                                 flowDescription, 
                                                                 simulationParameters)
{
  rnExact = PETSC_NULL;
} // ConvectiveTerm


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
PetscErrorCode ConvectiveTerm<dim>::initializeFluxes()
{
  return 0;
} // initializeFluxes


// two-dimensional specialization
template <>
PetscErrorCode ConvectiveTerm<2>::initializeFluxes()
{
  PetscErrorCode ierr;

  PetscInt i, j,           // loop indices
           m, n,           // local number of nodes along each direction
           mstart, nstart; // starting indices

  PetscReal x, y; // coordinates of a node

  Vec qxGlobal, qyGlobal;
  ierr = DMCompositeGetAccess(qPack, q, &qxGlobal, &qyGlobal); CHKERRQ(ierr);

  // fluxes in x-direction
  PetscReal **qx;
  ierr = DMDAVecGetArray(uda, qxGlobal, &qx); CHKERRQ(ierr);
  ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
  PetscReal u;  // value of u-velocity
  for (j=nstart; j<nstart+n; j++)
  {
    for (i=mstart; i<mstart+m; i++)
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
  ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
  PetscReal v;  // value of v-velocity
  for (j=nstart; j<nstart+n; j++)
  {
    for (i=mstart; i<mstart+m; i++)
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
} // initializeFluxes


// three-dimensional specialization
template <>
PetscErrorCode ConvectiveTerm<3>::initializeFluxes()
{
  PetscErrorCode ierr;

  PetscInt i, j, k,                // loop indices
           m, n, p,                // local number of nodes along each direction
           mstart, nstart, pstart; // starting indices

  PetscReal x, y, z; // coordinates of a node

  Vec qxGlobal, qyGlobal, qzGlobal;
  ierr = DMCompositeGetAccess(qPack, q, &qxGlobal, &qyGlobal, &qzGlobal); CHKERRQ(ierr);

  // fluxes in x-direction
  PetscReal ***qx;
  ierr = DMDAVecGetArray(uda, qxGlobal, &qx); CHKERRQ(ierr);
  ierr = DMDAGetCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  PetscReal u;  // value of u-velocity
  for (k=pstart; k<pstart+p; k++)
  {
    for (j=nstart; j<nstart+n; j++)
    {
      for (i=mstart; i<mstart+m; i++)
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
  ierr = DMDAGetCorners(vda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  PetscReal v;  // value of v-velocity
  for (k=pstart; k<pstart+p; k++)
  {
    for (j=nstart; j<nstart+n; j++)
    {
      for (i=mstart; i<mstart+m; i++)
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
  ierr = DMDAGetCorners(wda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  PetscReal w;  // value of w-velocity
  for (k=pstart; k<pstart+p; k++)
  {
    for (j=nstart; j<nstart+n; j++)
    {
      for (i=mstart; i<mstart+m; i++)
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
} // initializeFluxes


/**
 * \brief Computes the exact solution explicit convective terms.
 */
template <PetscInt dim>
PetscErrorCode ConvectiveTerm<dim>::calculateExactSolution()
{
  return 0;
} // calculateExactSolution


// two-dimensional specialization
template <>
PetscErrorCode ConvectiveTerm<2>::calculateExactSolution()
{
  PetscErrorCode ierr;
  
  PetscInt i, j,           // loop indices
           m, n,           // local number of nodes along each direction
           mstart, nstart; // starting indices
 
  PetscReal x, y,       // coordinates of a node
            velocity,   // velocity value at the node
            convection; // convective term at the node
  
  PetscReal dt = parameters->dt,                            // time-increment
            gamma = parameters->convection.coefficients[1]; // explicit time-scheme coefficient for convective terms

  ierr = VecDuplicate(q, &rnExact); CHKERRQ(ierr);
  Vec rnExactXGlobal, rnExactYGlobal;
  ierr = DMCompositeGetAccess(qPack, rnExact, &rnExactXGlobal, &rnExactYGlobal); CHKERRQ(ierr);

  // exact explicit terms for u-nodes
  PetscReal **rnExactX;
  ierr = DMDAVecGetArray(uda, rnExactXGlobal, &rnExactX);
  ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
  for (j=nstart; j<nstart+n; j++)
  {
    for (i=mstart; i<mstart+m; i++)
    {
      x = mesh->x[i+1];
      y = 0.5*(mesh->y[j]+mesh->y[j+1]);
      velocity = sin(PETSC_PI*x)*sin(PETSC_PI*y) + 1.0;
      convection = 2.0*PETSC_PI*(cos(PETSC_PI*x)*sin(PETSC_PI*y))*(sin(PETSC_PI*x)*sin(PETSC_PI*y)+1.0) \
                   + PETSC_PI*(sin(PETSC_PI*x)*cos(PETSC_PI*y))*(2.0*sin(PETSC_PI*x)*sin(PETSC_PI*y)+1.0);
      rnExactX[j][i] = velocity/dt - gamma*convection;
    }
  }
  ierr = DMDAVecRestoreArray(uda, rnExactXGlobal, &rnExactX); CHKERRQ(ierr);

  // exact explicit terms for v-nodes
  PetscReal **rnExactY;
  ierr = DMDAVecGetArray(vda, rnExactYGlobal, &rnExactY);
  ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
  for (j=nstart; j<nstart+n; j++)
  {
    for (i=mstart; i<mstart+m; i++)
    {
      x = 0.5*(mesh->x[i]+mesh->x[i+1]);
      y = mesh->y[j+1];
      velocity = sin(PETSC_PI*x)*sin(PETSC_PI*y);
      convection = PETSC_PI*(cos(PETSC_PI*x)*sin(PETSC_PI*y))*(2.0*sin(PETSC_PI*x)*sin(PETSC_PI*y)+1.0) \
                   + 2.0*PETSC_PI*(sin(PETSC_PI*x)*cos(PETSC_PI*y))*(sin(PETSC_PI*x)*sin(PETSC_PI*y));
      rnExactY[j][i] = velocity/dt - gamma*convection;
    }
  }
  ierr = DMDAVecRestoreArray(vda, rnExactYGlobal, &rnExactY); CHKERRQ(ierr);

  ierr = DMCompositeRestoreAccess(qPack, rnExact, &rnExactXGlobal, &rnExactYGlobal); CHKERRQ(ierr);

  return ierr;
} // calculateExactSolution


// three-dimensional specialization
template <>
PetscErrorCode ConvectiveTerm<3>::calculateExactSolution()
{
  PetscErrorCode ierr;
  
  PetscInt i, j, k,                // loop indices
           m, n, p,                // local number of nodes along each direction
           mstart, nstart, pstart; // starting indices
  
  PetscReal x, y, z,    // coordinates of a node
            velocity,   // velocity value at the node
            convection; // convective term at the node
  
  PetscReal dt = parameters->dt,                           // time-increment
            gamma = parameters->convection.coefficients[1]; // explicit time-scheme coefficient for convective terms

  ierr = VecDuplicate(q, &rnExact); CHKERRQ(ierr);
  Vec rnExactXGlobal, rnExactYGlobal, rnExactZGlobal;
  ierr = DMCompositeGetAccess(qPack, rnExact, &rnExactXGlobal, &rnExactYGlobal, &rnExactZGlobal); CHKERRQ(ierr);

  // exact explicit terms for u-nodes
  PetscReal ***rnExactX;
  ierr = DMDAVecGetArray(uda, rnExactXGlobal, &rnExactX);
  ierr = DMDAGetCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  for (k=pstart; k<pstart+p; k++)
  {
    for (j=nstart; j<nstart+n; j++)
    {
      for (i=mstart; i<mstart+m; i++)
      {
        x = mesh->x[i+1];
        y = 0.5*(mesh->y[j]+mesh->y[j+1]);
        z = 0.5*(mesh->z[k]+mesh->z[k+1]);
        velocity = sin(PETSC_PI*x)*sin(PETSC_PI*y)*sin(PETSC_PI*z) + 1.0;
        convection = 2.0*PETSC_PI*(cos(PETSC_PI*x)*sin(PETSC_PI*y)*sin(PETSC_PI*z))*(sin(PETSC_PI*x)*sin(PETSC_PI*y)*sin(PETSC_PI*z)+1.0) \
                     + PETSC_PI*(sin(PETSC_PI*x)*cos(PETSC_PI*y)*sin(PETSC_PI*z))*(2.0*sin(PETSC_PI*x)*sin(PETSC_PI*y)*sin(PETSC_PI*z)+1.0) \
                     + PETSC_PI*(sin(PETSC_PI*x)*sin(PETSC_PI*y)*cos(PETSC_PI*z))*(2.0*sin(PETSC_PI*x)*sin(PETSC_PI*y)*sin(PETSC_PI*z)+1.0);
        rnExactX[k][j][i] = velocity/dt - gamma*convection;
      }
    }
  }
  ierr = DMDAVecRestoreArray(uda, rnExactXGlobal, &rnExactX); CHKERRQ(ierr);

  // exact explicit terms for v-nodes
  PetscReal ***rnExactY;
  ierr = DMDAVecGetArray(vda, rnExactYGlobal, &rnExactY);
  ierr = DMDAGetCorners(vda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  for (k=pstart; k<pstart+p; k++)
  {
    for (j=nstart; j<nstart+n; j++)
    {
      for (i=mstart; i<mstart+m; i++)
      {
        x = 0.5*(mesh->x[i]+mesh->x[i+1]);
        y = mesh->y[j+1];
        z = 0.5*(mesh->z[k]+mesh->z[k+1]);
        velocity = sin(PETSC_PI*x)*sin(PETSC_PI*y)*sin(PETSC_PI*z);
        convection = PETSC_PI*(cos(PETSC_PI*x)*sin(PETSC_PI*y)*sin(PETSC_PI*z))*(2.0*sin(PETSC_PI*x)*sin(PETSC_PI*y)*sin(PETSC_PI*z)+1.0) \
                     + 2.0*PETSC_PI*(sin(PETSC_PI*x)*cos(PETSC_PI*y)*sin(PETSC_PI*z))*(sin(PETSC_PI*x)*sin(PETSC_PI*y)*sin(PETSC_PI*z)) \
                     + 2.0*PETSC_PI*(sin(PETSC_PI*x)*sin(PETSC_PI*y)*cos(PETSC_PI*z))*(sin(PETSC_PI*x)*sin(PETSC_PI*y)*sin(PETSC_PI*z));
        rnExactY[k][j][i] = velocity/dt - gamma*convection;
      }
    }
  }
  ierr = DMDAVecRestoreArray(vda, rnExactYGlobal, &rnExactY); CHKERRQ(ierr);

  // exact explicit terms for w-nodes
  PetscReal ***rnExactZ;
  ierr = DMDAVecGetArray(wda, rnExactZGlobal, &rnExactZ);
  ierr = DMDAGetCorners(wda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  for (k=pstart; k<pstart+p; k++)
  {
    for (j=nstart; j<nstart+n; j++)
    {
      for (i=mstart; i<mstart+m; i++)
      {
        x = 0.5*(mesh->x[i]+mesh->x[i+1]);
        y = 0.5*(mesh->y[j]+mesh->y[j+1]);
        z = mesh->z[k+1];
        velocity = sin(PETSC_PI*x)*sin(PETSC_PI*y)*sin(PETSC_PI*z);
        convection = PETSC_PI*(cos(PETSC_PI*x)*sin(PETSC_PI*y)*sin(PETSC_PI*z))*(2.0*sin(PETSC_PI*x)*sin(PETSC_PI*y)*sin(PETSC_PI*z)+1.0) \
                     + 2.0*PETSC_PI*(sin(PETSC_PI*x)*cos(PETSC_PI*y)*sin(PETSC_PI*z))*(sin(PETSC_PI*x)*sin(PETSC_PI*y)*sin(PETSC_PI*z)) \
                     + 2.0*PETSC_PI*(sin(PETSC_PI*x)*sin(PETSC_PI*y)*cos(PETSC_PI*z))*(sin(PETSC_PI*x)*sin(PETSC_PI*y)*sin(PETSC_PI*z));
        rnExactZ[k][j][i] = velocity/dt - gamma*convection;
      }
    }
  }
  ierr = DMDAVecRestoreArray(wda, rnExactZGlobal, &rnExactZ); CHKERRQ(ierr);

  ierr = DMCompositeRestoreAccess(qPack, rnExact, &rnExactXGlobal, &rnExactYGlobal, &rnExactZGlobal); CHKERRQ(ierr);

  return ierr;
} // calculateExactSolution


/**
 * \brief Computes the relative L2 norm of the difference between the numerical
 *        and the exact solution of the explicit terms.
 */
template<PetscInt dim>
PetscErrorCode ConvectiveTerm<dim>::calculateRelativeError()
{
  PetscErrorCode ierr;

  ierr = VecAXPY(NavierStokesSolver<dim>::rn, -1.0, rnExact);
  PetscReal l2NormDifference;
  ierr = VecNorm(NavierStokesSolver<dim>::rn, NORM_2, &l2NormDifference); CHKERRQ(ierr);
  PetscReal l2NormExact;
  ierr = VecNorm(rnExact, NORM_2, &l2NormExact); CHKERRQ(ierr);
  relativeError = l2NormDifference/l2NormExact;

  return ierr;
} // calculateRelativeError


/**
 * \brief Writes the number of cells and the relative error in a file.
 *
 * In the file already exists, the values are appended; otherwise the file 
 * is created.
 */
template <PetscInt dim>
PetscErrorCode ConvectiveTerm<dim>::writeRelativeError()
{
  PetscErrorCode ierr;

  PetscInt rank;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  if (rank == 0)
  {
    std::stringstream out;
    out << "./data/relativeErrors" << dim << "d.dat";
    std::ofstream fileStream;
    fileStream.open(out.str().c_str(), std::ios::out | std::ios::app);
    PetscInt nx = NavierStokesSolver<dim>::mesh->nx;
    PetscReal xStart = NavierStokesSolver<dim>::mesh->x[0],
              xEnd = NavierStokesSolver<dim>::mesh->x[nx];
    PetscReal h = (xEnd-xStart)/nx;
    fileStream << h << "\t" << relativeError << std::endl;
    fileStream.close();
  }

  return ierr;
} // writeRelativeError


/**
 * \brief Frees memory to avoid memory leaks
 */
template <PetscInt dim>
PetscErrorCode ConvectiveTerm<dim>::finalize()
{
  PetscErrorCode ierr;

  if (rnExact != PETSC_NULL)
  {
    ierr = VecDestroy(&rnExact); CHKERRQ(ierr);
  }
  ierr = NavierStokesSolver<dim>::finalize(); CHKERRQ(ierr);
  
  return ierr;
} // finalize


// template class specialization
template class ConvectiveTerm<2>;
template class ConvectiveTerm<3>;