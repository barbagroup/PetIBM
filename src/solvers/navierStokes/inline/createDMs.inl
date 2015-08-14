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
  DMBoundaryType dmBoundaryX = (flow->boundaries[XMINUS][U].type == PERIODIC) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED, 
                 dmBoundaryY = (flow->boundaries[YMINUS][U].type == PERIODIC) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;
    
  // create DMDA object for pressure
  numX = mesh->nx;
  numY = mesh->ny;
  ierr = DMDACreate2d(PETSC_COMM_WORLD, 
                      dmBoundaryX, dmBoundaryY, 
                      DMDA_STENCIL_STAR, 
                      numX, numY, PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL, 
                      &pda); CHKERRQ(ierr);

  // create DMDA objects for fluxes using the one for pressure
  const PetscInt *lxp, *lyp;
  ierr = DMDAGetOwnershipRanges(pda, &lxp, &lyp, NULL); CHKERRQ(ierr);
  PetscInt m, n;
  ierr = DMDAGetInfo(pda, NULL, NULL, NULL, NULL, &m, &n, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
  // fluxes in x-direction
  PetscInt *lxu, *lyu;
  ierr = PetscMalloc(m*sizeof(*lxu), &lxu); CHKERRQ(ierr);
  ierr = PetscMalloc(n*sizeof(*lyu), &lyu); CHKERRQ(ierr);
  ierr = PetscMemcpy(lxu, lxp, m*sizeof(*lxu)); CHKERRQ(ierr);
  ierr = PetscMemcpy(lyu, lyp, n*sizeof(*lyu)); CHKERRQ(ierr);
  numX = mesh->nx;
  numY = mesh->ny;
  if (flow->boundaries[XMINUS][U].type != PERIODIC)
  {
    lxu[m-1]--;
    numX--;
  }
  ierr = DMDACreate2d(PETSC_COMM_WORLD, 
                      dmBoundaryX, dmBoundaryY, 
                      DMDA_STENCIL_BOX, 
                      numX, numY, m, n, 1, 1, lxu, lyu, 
                      &uda); CHKERRQ(ierr);
  ierr = PetscFree(lxu); CHKERRQ(ierr);
  ierr = PetscFree(lyu); CHKERRQ(ierr);
  // fluxes in y-direction
  PetscInt *lxv, *lyv;
  ierr = PetscMalloc(m*sizeof(*lxv), &lxv); CHKERRQ(ierr);
  ierr = PetscMalloc(n*sizeof(*lyv), &lyv); CHKERRQ(ierr);
  ierr = PetscMemcpy(lxv, lxp, m*sizeof(*lxv)); CHKERRQ(ierr);
  ierr = PetscMemcpy(lyv, lyp, n*sizeof(*lyv)); CHKERRQ(ierr);
  numX = mesh->nx;
  numY = mesh->ny;
  if (flow->boundaries[YMINUS][V].type != PERIODIC)
  {
    lyv[n-1]--;
    numY--;
  }
  ierr = DMDACreate2d(PETSC_COMM_WORLD, 
                      dmBoundaryX, dmBoundaryY, 
                      DMDA_STENCIL_BOX, 
                      numX, numY, m, n, 1, 1, lxv, lyv, 
                      &vda); CHKERRQ(ierr);
  ierr = PetscFree(lxv); CHKERRQ(ierr);
  ierr = PetscFree(lyv); CHKERRQ(ierr);

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
  DMBoundaryType dmBoundaryX = (flow->boundaries[XMINUS][U].type == PERIODIC) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED, 
                 dmBoundaryY = (flow->boundaries[YMINUS][U].type == PERIODIC) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED, 
                 dmBoundaryZ = (flow->boundaries[ZMINUS][U].type == PERIODIC) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED;

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
  const PetscInt *lxp, *lyp, *lzp;
  ierr = DMDAGetOwnershipRanges(pda, &lxp, &lyp, &lzp); CHKERRQ(ierr);
  PetscInt m, n, p;
  ierr = DMDAGetInfo(pda, NULL, NULL, NULL, NULL, &m, &n, &p, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
  // fluxes in x-direction
  PetscInt *lxu, *lyu, *lzu;
  ierr = PetscMalloc(m*sizeof(*lxu), &lxu); CHKERRQ(ierr);
  ierr = PetscMalloc(n*sizeof(*lyu), &lyu); CHKERRQ(ierr);
  ierr = PetscMalloc(p*sizeof(*lzu), &lzu); CHKERRQ(ierr);
  ierr = PetscMemcpy(lxu, lxp, m*sizeof(*lxu)); CHKERRQ(ierr);
  ierr = PetscMemcpy(lyu, lyp, n*sizeof(*lyu)); CHKERRQ(ierr);
  ierr = PetscMemcpy(lzu, lzp, p*sizeof(*lzu)); CHKERRQ(ierr);
  numX = mesh->nx;
  numY = mesh->ny;
  numZ = mesh->nz;
  if (flow->boundaries[XMINUS][U].type != PERIODIC)
  {
    lxu[m-1]--;
    numX--;
  }
  ierr = DMDACreate3d(PETSC_COMM_WORLD, 
                      dmBoundaryX, dmBoundaryY, dmBoundaryZ, 
                      DMDA_STENCIL_BOX, 
                      numX, numY, numZ, m, n, p, 1, 1, lxu, lyu, lzu, 
                      &uda); CHKERRQ(ierr);
  ierr = PetscFree(lxu); CHKERRQ(ierr);
  ierr = PetscFree(lyu); CHKERRQ(ierr);
  ierr = PetscFree(lzu); CHKERRQ(ierr);
  // fluxes in y-direction
  PetscInt *lxv, *lyv, *lzv;
  ierr = PetscMalloc(m*sizeof(*lxv), &lxv); CHKERRQ(ierr);
  ierr = PetscMalloc(n*sizeof(*lyv), &lyv); CHKERRQ(ierr);
  ierr = PetscMalloc(p*sizeof(*lzv), &lzv); CHKERRQ(ierr);
  ierr = PetscMemcpy(lxv, lxp, m*sizeof(*lxv)); CHKERRQ(ierr);
  ierr = PetscMemcpy(lyv, lyp, n*sizeof(*lyv)); CHKERRQ(ierr);
  ierr = PetscMemcpy(lzv, lzp, p*sizeof(*lzv)); CHKERRQ(ierr);
  numX = mesh->nx;
  numY = mesh->ny;
  numZ = mesh->nz;
  if (flow->boundaries[YMINUS][V].type != PERIODIC)
  {
    lyv[n-1]--;
    numY--;
  }
  ierr = DMDACreate3d(PETSC_COMM_WORLD, 
                      dmBoundaryX, dmBoundaryY, dmBoundaryZ, 
                      DMDA_STENCIL_BOX, 
                      numX, numY, numZ, m, n, p, 1, 1, lxv, lyv, lzv, 
                      &vda); CHKERRQ(ierr);
  ierr = PetscFree(lxv); CHKERRQ(ierr);
  ierr = PetscFree(lyv); CHKERRQ(ierr);
  ierr = PetscFree(lzv); CHKERRQ(ierr);
  // fluxes in z-direction
  PetscInt *lxw, *lyw, *lzw;
  ierr = PetscMalloc(m*sizeof(*lxw), &lxw); CHKERRQ(ierr);
  ierr = PetscMalloc(n*sizeof(*lyw), &lyw); CHKERRQ(ierr);
  ierr = PetscMalloc(p*sizeof(*lzw), &lzw); CHKERRQ(ierr);
  ierr = PetscMemcpy(lxw, lxp, m*sizeof(*lxw)); CHKERRQ(ierr);
  ierr = PetscMemcpy(lyw, lyp, n*sizeof(*lyw)); CHKERRQ(ierr);
  ierr = PetscMemcpy(lzw, lzp, p*sizeof(*lzw)); CHKERRQ(ierr);
  numX = mesh->nx;
  numY = mesh->ny;
  numZ = mesh->nz;
  if (flow->boundaries[ZMINUS][W].type != PERIODIC)
  {
    lzw[p-1]--;
    numZ--;
  }
  ierr = DMDACreate3d(PETSC_COMM_WORLD, 
                      dmBoundaryX, dmBoundaryY, dmBoundaryZ, 
                      DMDA_STENCIL_BOX, 
                      numX, numY, numZ, m, n, p, 1, 1, lxw, lyw, lzw, 
                      &wda); CHKERRQ(ierr);
  ierr = PetscFree(lxw); CHKERRQ(ierr);
  ierr = PetscFree(lyw); CHKERRQ(ierr);
  ierr = PetscFree(lzw); CHKERRQ(ierr);

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