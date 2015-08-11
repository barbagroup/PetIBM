/***************************************************************************//**
 * \file initializeMeshSpacings.inl
 * \author Anush Krishnan (anush@bu.edu), Olivier Mesnard (mesnardo@gwu.edu)
 * \brief Implementation of the method `initializeMeshSpacings`.
 */


/**
 * \brief Computes the distance between consecutive fluxes in each direction
 *        and along each direction.
 *
 * The cell widths stored in the CartesianMesh object `mesh` refer to the 
 * widths of the cells of the grid used to discretize the domain. But we require
 * the distance between consecutive velocity flux locations for certain 
 * functions (e.g. when we assemble the matrix `A`, or update the ghost cells on 
 * the domain boundary.) This function calculates the spacings between the 
 * points where the velocity fluxes are computed. 
 *
 * At locations near the domain boundaries, the distance from the velocity flux
 * to the boundary is calculated. In the case of periodic domains, the distance 
 * between the velocity fluxes at the opposite edges are calculated, assuming
 * that the domain has been wrapped around.
 */
template <PetscInt dim>
void NavierStokesSolver<dim>::initializeMeshSpacings()
{
  PetscInt i, j, k; // loop indices
  PetscInt numX, // number of values along x-direction
           numY, // number of values along y-direction
           numZ; // number of values along z-direction

  // look for periodic boundary conditions in x- and y- directions
  bool periodicX = (flow->boundaries[XMINUS][0].type == PERIODIC);
  bool periodicY = (flow->boundaries[YMINUS][0].type == PERIODIC); 

  // spacings along x-direction between fluxes in x-direction
  numX = (periodicX) ? mesh->nx : mesh->nx-1;
  dxU.resize(numX+1);
  for (i=0; i<numX; i++)
  {
    dxU[i] = mesh->dx[i];
  }
  dxU[numX] = mesh->dx[0];
  // spacings along y-direction between fluxes in x-direction
  numY = mesh->ny;
  dyU.resize(numY+1);
  dyU[0] = (periodicY) ? 0.5*(mesh->dy[mesh->ny-1]+mesh->dy[0]) : 0.5*(mesh->dy[0]);
  for (j=1; j<numY; j++)
  {
    dyU[j] = 0.5 * (mesh->dy[j-1] + mesh->dy[j]);
  }
  dyU[numY] = (periodicY) ? 0.5*(mesh->dy[mesh->ny-1]+mesh->dy[0]) : 0.5*(mesh->dy[mesh->ny-1]);

  // spacings along x-direction for fluxes in y-direction
  numX = mesh->nx;
  dxV.resize(numX+1);
  dxV[0] = (periodicX) ? 0.5*(mesh->dx[mesh->nx-1]+mesh->dx[0]) : 0.5*(mesh->dx[0]);
  for (i=1; i<numX; i++)
  {
    dxV[i] = 0.5 * (mesh->dx[i-1] + mesh->dx[i]);
  }
  dxV[numX] = (periodicX) ? 0.5*(mesh->dx[mesh->nx-1]+mesh->dx[0]) : 0.5*(mesh->dx[mesh->nx-1]);
  // spacings along y-direction for fluxes in y-direction
  numY = (periodicY) ? mesh->ny : mesh->ny-1;
  dyV.resize(numY+1);
  for (j=0; j<numY; j++)
  {
    dyV[j] = mesh->dy[j];
  }
  dyV[numY] = mesh->dy[0];

  if (dim == 3)
  {
    // look for periodic boundary conditions in z-direction
    bool periodicZ = (flow->boundaries[ZMINUS][0].type == PERIODIC);

    // spacings along z-direction for fluxes in x- and y- directions
    numZ = mesh->nz;
    dzV.resize(numZ+1);
    dzU[0] = (periodicZ) ? 0.5*(mesh->dz[mesh->nz-1]+mesh->dz[0]) : 0.5*(mesh->dz[0]);
    dzV[0] = dzU[0];
    for (k=1; k<numZ; k++)
    {
      dzU[k] = 0.5 * (mesh->dz[k-1] + mesh->dz[k]);
      dzV[k] = dzU[k];
    }
    dzU[numZ] = (periodicZ) ? 0.5*(mesh->dz[mesh->nz-1]+mesh->dz[0]) : 0.5*(mesh->dz[mesh->nz-1]);
    dzV[numZ] = dzU[numZ];

    // spacings along x-direction for fluxes in z-direction
    for (i=0; i<dxV.size(); i++)
    {
      dxW[i] = dxV[i];
    }
    // spacings along y-direction for fluxes in z-direction
    for (j=0; j<dyU.size(); j++)
    {
      dyW[j] = dyU[j];
    }
    // spacings along z-direction for fluxes in z-direction
    numZ = (periodicZ) ? mesh->nz : mesh->nz-1;
    dzW.resize(numZ+1);
    for (k=0; k<numZ; k++)
    {
      dzW[k] = mesh->dz[k];
    }
    dzV[numZ] = mesh->dz[0];
  }
  
} // initializeMeshSpacings