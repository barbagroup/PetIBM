/***************************************************************************//**
* Create the data structures used to store information about the distributed
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
}

template <>
PetscErrorCode NavierStokesSolver<2>::createDMs()
{
	PetscErrorCode    ierr;
	PetscInt          m, n;
	const PetscInt    *lxp, *lyp;
	PetscInt          *lxu, *lyu, *lxv, *lyv;
	PetscInt          numX, numY;
	DMDABoundaryType  bx, by;
	
	// set boundary types
	bx = (flowDesc->bc[0][XPLUS].type == PERIODIC)? DMDA_BOUNDARY_PERIODIC : DMDA_BOUNDARY_GHOSTED;
	by = (flowDesc->bc[0][YPLUS].type == PERIODIC)? DMDA_BOUNDARY_PERIODIC : DMDA_BOUNDARY_GHOSTED;	
	
	// Create distributed array data structures
	// pressure
	numX = mesh->nx;
	numY = mesh->ny;
	ierr = DMDACreate2d(PETSC_COMM_WORLD, bx, by, DMDA_STENCIL_STAR, numX, numY, PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL, &pda); CHKERRQ(ierr);
	ierr = DMDAGetOwnershipRanges(pda, &lxp, &lyp, NULL); CHKERRQ(ierr);
	ierr = DMDAGetInfo(pda, NULL, NULL, NULL, NULL, &m, &n, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
	
	// packed DMs
	ierr = DMCompositeCreate(PETSC_COMM_WORLD, &qPack); CHKERRQ(ierr);
	ierr = DMCompositeCreate(PETSC_COMM_WORLD, &lambdaPack); CHKERRQ(ierr);
	ierr = DMCompositeAddDM(lambdaPack, pda); CHKERRQ(ierr);
	
	// x-velocity
	ierr = PetscMalloc(m*sizeof(*lxu), &lxu); CHKERRQ(ierr);
	ierr = PetscMalloc(n*sizeof(*lyu), &lyu); CHKERRQ(ierr);
	ierr = PetscMemcpy(lxu, lxp, m*sizeof(*lxu)); CHKERRQ(ierr);
	ierr = PetscMemcpy(lyu, lyp, n*sizeof(*lyu)); CHKERRQ(ierr);
	numX = mesh->nx;
	numY = mesh->ny;
	if(flowDesc->bc[0][XPLUS].type != PERIODIC)
	{
		lxu[m-1]--;
		numX = mesh->nx-1;
	}
	ierr = DMDACreate2d(PETSC_COMM_WORLD, bx, by, DMDA_STENCIL_BOX, numX, numY, m, n, 1, 1, lxu, lyu, &uda); CHKERRQ(ierr);
	ierr = DMCompositeAddDM(qPack, uda); CHKERRQ(ierr);
	
	// y-velocity
	ierr = PetscMalloc(m*sizeof(*lxv), &lxv); CHKERRQ(ierr);
	ierr = PetscMalloc(n*sizeof(*lyv), &lyv); CHKERRQ(ierr);
	ierr = PetscMemcpy(lxv, lxp, m*sizeof(*lxv)); CHKERRQ(ierr);
	ierr = PetscMemcpy(lyv, lyp, n*sizeof(*lyv)); CHKERRQ(ierr);
	numX = mesh->nx;
	numY = mesh->ny;
	if(flowDesc->bc[1][YPLUS].type != PERIODIC)
	{
		lyv[n-1]--;
		numY = mesh->ny-1;
	}
	ierr = DMDACreate2d(PETSC_COMM_WORLD, bx, by, DMDA_STENCIL_BOX, numX, numY, m, n, 1, 1, lxv, lyv, &vda); CHKERRQ(ierr);
	ierr = DMCompositeAddDM(qPack, vda); CHKERRQ(ierr);
	
	PetscFree(lxu);
	PetscFree(lyu);	
	PetscFree(lxv);
	PetscFree(lyv);

	return 0;
}

template <>
PetscErrorCode NavierStokesSolver<3>::createDMs()
{
	PetscErrorCode    ierr;
	PetscInt          m, n, p;
	const PetscInt    *lxp, *lyp, *lzp;
	PetscInt          *lxu, *lyu, *lzu;
	PetscInt          *lxv, *lyv, *lzv;
	PetscInt          *lxw, *lyw, *lzw;
	PetscInt          numX, numY, numZ;
	DMDABoundaryType  bx, by, bz;
	
	// set boundary types
	bx = (flowDesc->bc[0][XPLUS].type == PERIODIC)? DMDA_BOUNDARY_PERIODIC : DMDA_BOUNDARY_GHOSTED;
	by = (flowDesc->bc[0][YPLUS].type == PERIODIC)? DMDA_BOUNDARY_PERIODIC : DMDA_BOUNDARY_GHOSTED;
	bz = (flowDesc->bc[0][ZPLUS].type == PERIODIC)? DMDA_BOUNDARY_PERIODIC : DMDA_BOUNDARY_GHOSTED;
	
	// Create distributed array data structures
	// pressure
	numX = mesh->nx;
	numY = mesh->ny;
	numZ = mesh->nz;
	ierr = DMDACreate3d(PETSC_COMM_WORLD, bx, by, bz, DMDA_STENCIL_STAR, numX, numY, numZ, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL, NULL, &pda); CHKERRQ(ierr);
	ierr = DMDAGetOwnershipRanges(pda, &lxp, &lyp, &lzp); CHKERRQ(ierr);
	ierr = DMDAGetInfo(pda, NULL, NULL, NULL, NULL, &m, &n, &p, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
	
	// packed DMs
	ierr = DMCompositeCreate(PETSC_COMM_WORLD, &qPack); CHKERRQ(ierr);
	ierr = DMCompositeCreate(PETSC_COMM_WORLD, &lambdaPack); CHKERRQ(ierr);
	ierr = DMCompositeAddDM(lambdaPack, pda); CHKERRQ(ierr);
	
	// x-velocity
	ierr = PetscMalloc(m*sizeof(*lxu), &lxu); CHKERRQ(ierr);
	ierr = PetscMalloc(n*sizeof(*lyu), &lyu); CHKERRQ(ierr);
	ierr = PetscMalloc(p*sizeof(*lzu), &lzu); CHKERRQ(ierr);
	ierr = PetscMemcpy(lxu, lxp, m*sizeof(*lxu)); CHKERRQ(ierr);
	ierr = PetscMemcpy(lyu, lyp, n*sizeof(*lyu)); CHKERRQ(ierr);
	ierr = PetscMemcpy(lzu, lzp, p*sizeof(*lzu)); CHKERRQ(ierr);
	numX = mesh->nx;
	numY = mesh->ny;
	numZ = mesh->nz;
	if(flowDesc->bc[0][XPLUS].type != PERIODIC)
	{
		lxu[m-1]--;
		numX = mesh->nx-1;
	}
	ierr = DMDACreate3d(PETSC_COMM_WORLD, bx, by, bz, DMDA_STENCIL_BOX, numX, numY, numZ, m, n, p, 1, 1, lxu, lyu, lzu, &uda); CHKERRQ(ierr);	
	ierr = DMCompositeAddDM(qPack, uda); CHKERRQ(ierr);
	
	// y-velocity
	ierr = PetscMalloc(m*sizeof(*lxv), &lxv); CHKERRQ(ierr);
	ierr = PetscMalloc(n*sizeof(*lyv), &lyv); CHKERRQ(ierr);
	ierr = PetscMalloc(p*sizeof(*lzv), &lzv); CHKERRQ(ierr);
	ierr = PetscMemcpy(lxv, lxp, m*sizeof(*lxv)); CHKERRQ(ierr);
	ierr = PetscMemcpy(lyv, lyp, n*sizeof(*lyv)); CHKERRQ(ierr);
	ierr = PetscMemcpy(lzv, lzp, p*sizeof(*lzv)); CHKERRQ(ierr);
	numX = mesh->nx;
	numY = mesh->ny;
	numZ = mesh->nz;
	if(flowDesc->bc[1][YPLUS].type != PERIODIC)
	{
		lyv[n-1]--;
		numY = mesh->ny-1;
	}
	ierr = DMDACreate3d(PETSC_COMM_WORLD, bx, by, bz, DMDA_STENCIL_BOX, numX, numY, numZ, m, n, p, 1, 1, lxv, lyv, lzv, &vda); CHKERRQ(ierr);
	ierr = DMCompositeAddDM(qPack, vda); CHKERRQ(ierr);
	
	// z-velocity
	ierr = PetscMalloc(m*sizeof(*lxw), &lxw); CHKERRQ(ierr);
	ierr = PetscMalloc(n*sizeof(*lyw), &lyw); CHKERRQ(ierr);
	ierr = PetscMalloc(p*sizeof(*lzw), &lzw); CHKERRQ(ierr);
	ierr = PetscMemcpy(lxw, lxp, m*sizeof(*lxw)); CHKERRQ(ierr);
	ierr = PetscMemcpy(lyw, lyp, n*sizeof(*lyw)); CHKERRQ(ierr);
	ierr = PetscMemcpy(lzw, lzp, p*sizeof(*lzw)); CHKERRQ(ierr);
	numX = mesh->nx;
	numY = mesh->ny;
	numZ = mesh->nz;
	if(flowDesc->bc[2][ZPLUS].type != PERIODIC)
	{
		lzw[p-1]--;
		numZ = mesh->nz-1;
	}
	ierr = DMDACreate3d(PETSC_COMM_WORLD, bx, by, bz, DMDA_STENCIL_BOX, numX, numY, numZ, m, n, p, 1, 1, lxw, lyw, lzw, &wda); CHKERRQ(ierr);
	ierr = DMCompositeAddDM(qPack, wda); CHKERRQ(ierr);
	
	PetscFree(lxu);
	PetscFree(lyu);
	PetscFree(lzu);
	PetscFree(lxv);
	PetscFree(lyv);
	PetscFree(lzv);
	PetscFree(lxw);
	PetscFree(lyw);
	PetscFree(lzw);

	return 0;
}
