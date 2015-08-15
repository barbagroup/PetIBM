/***************************************************************************//**
 * \file createDMs.inl
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the method `createDMs` of the class `NavierStokesSolver`.
 */


/**
 * \brief Creates the data structures used to store information about the distributed 
 * arrays used to store the velocity fluxes and the pressure.
 *
 * These flow variables are stored on a staggered grid, with the velocity fluxes 
 * stored at the faces of the cells, and the pressures at the centers. Hence, 
 * the sizes of the arrays of each of the flow variables are:
 * \f[ u : (nx-1) \times ny \times nz \f]
 * \f[ v : nx \times (ny-1) \times nz \f]
 * \f[ w : nx \times ny \times (nz-1) \f]
 * \f[ p : nx \times ny \times nz \f]
 *
 * In periodic domains, an extra layer of flux variables is required to store 
 * the values on the periodic boundary. This is taken into account when 
 * creating the distributed arrays.
 *
 * The vector used to store the velocity fluxes is a composite of the individual 
 * vectors that store the fluxes in each cartesian direction.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::createDMs()
{
  return 0;
} // createDMs


// two-dimensional specialization
template <>
PetscErrorCode NavierStokesSolver<2>::createDMs()
{
  PetscErrorCode ierr;
  
  PetscInt numX, numY;
    
  // set boundary types (periodic or ghosted)
  DMBoundaryType dmBoundaryX = (flow->boundaries[XMINUS][0].type == PERIODIC) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED, 
                 dmBoundaryY = (flow->boundaries[YMINUS][0].type == PERIODIC) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
    
  // create DMDA object for pressure
  numX = mesh->nx;
  numY = mesh->ny;
  ierr = DMDACreate2d(PETSC_COMM_WORLD, 
                      dmBoundaryX, dmBoundaryY, 
                      DMDA_STENCIL_STAR, 
                      numX, numY, PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL, 
                      &pda); CHKERRQ(ierr);

  // create DMDA objects for fluxes using the one for pressure
  const PetscInt *plx, *ply;
  ierr = DMDAGetOwnershipRanges(pda, &plx, &ply, NULL); CHKERRQ(ierr);
  PetscInt m, n;
  ierr = DMDAGetInfo(pda, NULL, NULL, NULL, NULL, &m, &n, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
  // fluxes in x-direction
  PetscInt *ulx, *uly;
  ierr = PetscMalloc(m*sizeof(*ulx), &ulx); CHKERRQ(ierr);
  ierr = PetscMalloc(n*sizeof(*uly), &uly); CHKERRQ(ierr);
  ierr = PetscMemcpy(ulx, plx, m*sizeof(*ulx)); CHKERRQ(ierr);
  ierr = PetscMemcpy(uly, ply, n*sizeof(*uly)); CHKERRQ(ierr);
  numX = mesh->nx;
  numY = mesh->ny;
  if (flow->boundaries[XMINUS][0].type != PERIODIC)
  {
    ulx[m-1]--;
    numX--;
  }
  ierr = DMDACreate2d(PETSC_COMM_WORLD, 
                      dmBoundaryX, dmBoundaryY, 
                      DMDA_STENCIL_BOX, 
                      numX, numY, m, n, 1, 1, ulx, uly, 
                      &uda); CHKERRQ(ierr);
  ierr = PetscFree(ulx); CHKERRQ(ierr);
  ierr = PetscFree(uly); CHKERRQ(ierr);
  // fluxes in y-direction
  PetscInt *vlx, *vly;
  ierr = PetscMalloc(m*sizeof(*vlx), &vlx); CHKERRQ(ierr);
  ierr = PetscMalloc(n*sizeof(*vly), &vly); CHKERRQ(ierr);
  ierr = PetscMemcpy(vlx, plx, m*sizeof(*vlx)); CHKERRQ(ierr);
  ierr = PetscMemcpy(vly, ply, n*sizeof(*vly)); CHKERRQ(ierr);
  numX = mesh->nx;
  numY = mesh->ny;
  if (flow->boundaries[YMINUS][1].type != PERIODIC)
  {
    vly[n-1]--;
    numY--;
  }
  ierr = DMDACreate2d(PETSC_COMM_WORLD, 
                      dmBoundaryX, dmBoundaryY, 
                      DMDA_STENCIL_BOX, 
                      numX, numY, m, n, 1, 1, vlx, vly, 
                      &vda); CHKERRQ(ierr);
  ierr = PetscFree(vlx); CHKERRQ(ierr);
  ierr = PetscFree(vly); CHKERRQ(ierr);

  // create lambda packer
  ierr = DMCompositeCreate(PETSC_COMM_WORLD, &lambdaPack); CHKERRQ(ierr);
  ierr = DMCompositeAddDM(lambdaPack, pda); CHKERRQ(ierr);

  // create fluxes packer
  ierr = DMCompositeCreate(PETSC_COMM_WORLD, &qPack); CHKERRQ(ierr);
  ierr = DMCompositeAddDM(qPack, uda); CHKERRQ(ierr);
  ierr = DMCompositeAddDM(qPack, vda); CHKERRQ(ierr);

  return 0;
} // createDMs


// three-dimensional specialization
template <>
PetscErrorCode NavierStokesSolver<3>::createDMs()
{
  PetscErrorCode ierr;

  PetscInt numX, numY, numZ;
  
  // set boundary types (periodic or ghosted)
  DMBoundaryType dmBoundaryX = (flow->boundaries[XMINUS][0].type == PERIODIC) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED, 
                 dmBoundaryY = (flow->boundaries[YMINUS][0].type == PERIODIC) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED, 
                 dmBoundaryZ = (flow->boundaries[ZMINUS][0].type == PERIODIC) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;

  // create DMDA object for pressure
  numX = mesh->nx;
  numY = mesh->ny;
  numZ = mesh->nz;
  ierr = DMDACreate3d(PETSC_COMM_WORLD, 
                      dmBoundaryX, dmBoundaryY, dmBoundaryZ, 
                      DMDA_STENCIL_STAR, 
                      numX, numY, numZ, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL, NULL, 
                      &pda); CHKERRQ(ierr);
  
  // create DMDA objects for fluxes from DMDA object for pressure
  const PetscInt *plx, *ply, *plz;
  ierr = DMDAGetOwnershipRanges(pda, &plx, &ply, &plz); CHKERRQ(ierr);
  PetscInt m, n, p;
  ierr = DMDAGetInfo(pda, NULL, NULL, NULL, NULL, &m, &n, &p, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
  // fluxes in x-direction
  PetscInt *ulx, *uly, *ulz;
  ierr = PetscMalloc(m*sizeof(*ulx), &ulx); CHKERRQ(ierr);
  ierr = PetscMalloc(n*sizeof(*uly), &uly); CHKERRQ(ierr);
  ierr = PetscMalloc(p*sizeof(*ulz), &ulz); CHKERRQ(ierr);
  ierr = PetscMemcpy(ulx, plx, m*sizeof(*ulx)); CHKERRQ(ierr);
  ierr = PetscMemcpy(uly, ply, n*sizeof(*uly)); CHKERRQ(ierr);
  ierr = PetscMemcpy(ulz, plz, p*sizeof(*ulz)); CHKERRQ(ierr);
  numX = mesh->nx;
  numY = mesh->ny;
  numZ = mesh->nz;
  if (flow->boundaries[XMINUS][0].type != PERIODIC)
  {
    ulx[m-1]--;
    numX--;
  }
  ierr = DMDACreate3d(PETSC_COMM_WORLD, 
                      dmBoundaryX, dmBoundaryY, dmBoundaryZ, 
                      DMDA_STENCIL_BOX, 
                      numX, numY, numZ, m, n, p, 1, 1, ulx, uly, ulz, 
                      &uda); CHKERRQ(ierr);
  ierr = PetscFree(ulx); CHKERRQ(ierr);
  ierr = PetscFree(uly); CHKERRQ(ierr);
  ierr = PetscFree(ulz); CHKERRQ(ierr);
  // fluxes in y-direction
  PetscInt *vlx, *vly, *vlz;
  ierr = PetscMalloc(m*sizeof(*vlx), &vlx); CHKERRQ(ierr);
  ierr = PetscMalloc(n*sizeof(*vly), &vly); CHKERRQ(ierr);
  ierr = PetscMalloc(p*sizeof(*vlz), &vlz); CHKERRQ(ierr);
  ierr = PetscMemcpy(vlx, plx, m*sizeof(*vlx)); CHKERRQ(ierr);
  ierr = PetscMemcpy(vly, ply, n*sizeof(*vly)); CHKERRQ(ierr);
  ierr = PetscMemcpy(vlz, plz, p*sizeof(*vlz)); CHKERRQ(ierr);
  numX = mesh->nx;
  numY = mesh->ny;
  numZ = mesh->nz;
  if (flow->boundaries[YMINUS][1].type != PERIODIC)
  {
    vly[n-1]--;
    numY--;
  }
  ierr = DMDACreate3d(PETSC_COMM_WORLD, 
                      dmBoundaryX, dmBoundaryY, dmBoundaryZ, 
                      DMDA_STENCIL_BOX, 
                      numX, numY, numZ, m, n, p, 1, 1, vlx, vly, vlz, 
                      &vda); CHKERRQ(ierr);
  ierr = PetscFree(vlx); CHKERRQ(ierr);
  ierr = PetscFree(vly); CHKERRQ(ierr);
  ierr = PetscFree(vlz); CHKERRQ(ierr);
  // fluxes in z-direction
  PetscInt *wlx, *wly, *wlz;
  ierr = PetscMalloc(m*sizeof(*wlx), &wlx); CHKERRQ(ierr);
  ierr = PetscMalloc(n*sizeof(*wly), &wly); CHKERRQ(ierr);
  ierr = PetscMalloc(p*sizeof(*wlz), &wlz); CHKERRQ(ierr);
  ierr = PetscMemcpy(wlx, plx, m*sizeof(*wlx)); CHKERRQ(ierr);
  ierr = PetscMemcpy(wly, ply, n*sizeof(*wly)); CHKERRQ(ierr);
  ierr = PetscMemcpy(wlz, plz, p*sizeof(*wlz)); CHKERRQ(ierr);
  numX = mesh->nx;
  numY = mesh->ny;
  numZ = mesh->nz;
  if (flow->boundaries[ZMINUS][2].type != PERIODIC)
  {
    wlz[p-1]--;
    numZ--;
  }
  ierr = DMDACreate3d(PETSC_COMM_WORLD, 
                      dmBoundaryX, dmBoundaryY, dmBoundaryZ, 
                      DMDA_STENCIL_BOX, 
                      numX, numY, numZ, m, n, p, 1, 1, wlx, wly, wlz, 
                      &wda); CHKERRQ(ierr);
  ierr = PetscFree(wlx); CHKERRQ(ierr);
  ierr = PetscFree(wly); CHKERRQ(ierr);
  ierr = PetscFree(wlz); CHKERRQ(ierr);

  // create lambda packer and add DMDA to it
  ierr = DMCompositeCreate(PETSC_COMM_WORLD, &lambdaPack); CHKERRQ(ierr);
  ierr = DMCompositeAddDM(lambdaPack, pda); CHKERRQ(ierr);

  // create fluxes packer and add DMDA to it
  ierr = DMCompositeCreate(PETSC_COMM_WORLD, &qPack); CHKERRQ(ierr);
  ierr = DMCompositeAddDM(qPack, uda); CHKERRQ(ierr);
  ierr = DMCompositeAddDM(qPack, vda); CHKERRQ(ierr);
  ierr = DMCompositeAddDM(qPack, wda); CHKERRQ(ierr);

  return 0;
} // createDMs