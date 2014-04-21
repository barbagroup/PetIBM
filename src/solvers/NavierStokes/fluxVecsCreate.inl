#include <petscdmcomposite.h>

template <>
void NavierStokesSolver<2>::fluxVecsCreate()
{
	PetscErrorCode    ierr;
	PetscInt          m, n;
	const PetscInt    *lxu, *lyu;
	PetscInt          *lxv, *lyv;
	PetscInt          numX, numY;
	DMDABoundaryType  bx, by;
	
	// set boundary types for velocity flux variables
	bx = (flowDesc->bc[0][XPLUS].type == PERIODIC)? DMDA_BOUNDARY_PERIODIC : DMDA_BOUNDARY_GHOSTED;
	by = (flowDesc->bc[0][YPLUS].type == PERIODIC)? DMDA_BOUNDARY_PERIODIC : DMDA_BOUNDARY_GHOSTED;
	
	// Create distributed array data structures
	// x-velocity
	numX = (flowDesc->bc[0][XPLUS].type == PERIODIC)? mesh->nx : mesh->nx-1;
	numY = mesh->ny;
	ierr = DMDACreate2d(PETSC_COMM_WORLD, bx, by, DMDA_STENCIL_BOX, numX, numY, PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL, &(uda)); CHKERRV(ierr);
	ierr = DMCompositeCreate(PETSC_COMM_WORLD, &(pack)); CHKERRV(ierr);
	ierr = DMCompositeAddDM(pack, uda); CHKERRV(ierr);
	ierr = DMDAGetOwnershipRanges(uda, &lxu, &lyu, NULL); CHKERRV(ierr);
	ierr = DMDAGetInfo(uda, NULL, NULL, NULL, NULL, &m, &n, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRV(ierr);
	
	// y-velocity
	ierr = PetscMalloc(m*sizeof(*lxv), &lxv); CHKERRV(ierr);
	ierr = PetscMalloc(n*sizeof(*lyv), &lyv); CHKERRV(ierr);
	ierr = PetscMemcpy(lxv, lxu, m*sizeof(*lxv)); CHKERRV(ierr);
	ierr = PetscMemcpy(lyv, lyu, n*sizeof(*lyv)); CHKERRV(ierr);
	if(flowDesc->bc[0][YPLUS].type != PERIODIC)
		lyv[n-1]--;
	if(numX == mesh->nx-1)
		lxv[m-1]++;
	numX = mesh->nx;
	numY = (flowDesc->bc[0][YPLUS].type == PERIODIC)? mesh->ny : mesh->ny-1;
	ierr = DMDACreate2d(PETSC_COMM_WORLD, bx, by, DMDA_STENCIL_BOX, numX, numY, m, n, 1, 1, lxv, lyv, &(vda)); CHKERRV(ierr);
	ierr = DMCompositeAddDM(pack, vda); CHKERRV(ierr);
		
	PetscFree(lxv);
	PetscFree(lyv);
	
	// create vectors
	// velocity fluxes
	ierr = DMCreateLocalVector(uda, &qxLocal); CHKERRV(ierr);
	ierr = DMCreateLocalVector(vda, &qyLocal); CHKERRV(ierr);
	
	// convection terms
	//ierr = DMCreateGlobalVector(pack, &rn); CHKERRV(ierr);
	//ierr = DMCreateGlobalVector(pack, &H); CHKERRV(ierr);
}

template <>
void NavierStokesSolver<3>::fluxVecsCreate()
{
	PetscErrorCode    ierr;
	PetscInt          m, n, p;
	const PetscInt    *lxu, *lyu, *lzu;
	PetscInt          *lxv, *lyv, *lzv;
	PetscInt          *lxw, *lyw, *lzw;
	PetscInt          numX, numY, numZ;
	DMDABoundaryType  bx, by, bz;
	
	// set boundary types for velocity variables
	bx = (flowDesc->bc[0][XPLUS].type == PERIODIC)? DMDA_BOUNDARY_PERIODIC : DMDA_BOUNDARY_GHOSTED;
	by = (flowDesc->bc[0][YPLUS].type == PERIODIC)? DMDA_BOUNDARY_PERIODIC : DMDA_BOUNDARY_GHOSTED;
	bz = (flowDesc->bc[0][ZPLUS].type == PERIODIC)? DMDA_BOUNDARY_PERIODIC : DMDA_BOUNDARY_GHOSTED;
	
	// Create distributed array data structures
	// x-velocity
	numX = (flowDesc->bc[0][XPLUS].type == PERIODIC)? mesh->nx : mesh->nx-1;
	numY = mesh->ny;
	numZ = mesh->nz;
	ierr = DMDACreate3d(PETSC_COMM_WORLD, bx, by, bz, DMDA_STENCIL_BOX, numX, numY, numZ, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL, NULL, &(uda)); CHKERRV(ierr);
	ierr = DMCompositeCreate(PETSC_COMM_WORLD, &(pack)); CHKERRV(ierr);
	ierr = DMCompositeAddDM(pack, uda); CHKERRV(ierr);
	ierr = DMDAGetOwnershipRanges(uda, &lxu, &lyu, &lzu); CHKERRV(ierr);
	ierr = DMDAGetInfo(uda, NULL, NULL, NULL, NULL, &m, &n, &p, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRV(ierr);
	
	// y-velocity
	ierr = PetscMalloc(m*sizeof(*lxv), &lxv); CHKERRV(ierr);
	ierr = PetscMalloc(n*sizeof(*lyv), &lyv); CHKERRV(ierr);
	ierr = PetscMalloc(p*sizeof(*lzv), &lzv); CHKERRV(ierr);
	ierr = PetscMemcpy(lxv, lxu, m*sizeof(*lxv)); CHKERRV(ierr);
	ierr = PetscMemcpy(lyv, lyu, n*sizeof(*lyv)); CHKERRV(ierr);
	ierr = PetscMemcpy(lzv, lzu, p*sizeof(*lzv)); CHKERRV(ierr);
	if(flowDesc->bc[0][YPLUS].type != PERIODIC)
		lyv[n-1]--;
	if(numX == mesh->nx-1)
		lxv[m-1]++;
	numX = mesh->nx;
	numY = (flowDesc->bc[0][YPLUS].type == PERIODIC)? mesh->ny : mesh->ny-1;
	numZ = mesh->nz;
	ierr = DMDACreate3d(PETSC_COMM_WORLD, bx, by, bz, DMDA_STENCIL_BOX, numX, numY, numZ, m, n, p, 1, 1, lxv, lyv, lzv, &(vda)); CHKERRV(ierr);
	ierr = DMCompositeAddDM(pack, vda); CHKERRV(ierr);
	
	// z-velocity
	ierr = PetscMalloc(m*sizeof(*lxw), &lxw); CHKERRV(ierr);
	ierr = PetscMalloc(n*sizeof(*lyw), &lyw); CHKERRV(ierr);
	ierr = PetscMalloc(p*sizeof(*lzw), &lzw); CHKERRV(ierr);
	ierr = PetscMemcpy(lxw, lxv, m*sizeof(*lxw)); CHKERRV(ierr);
	ierr = PetscMemcpy(lyw, lyv, n*sizeof(*lyw)); CHKERRV(ierr);
	ierr = PetscMemcpy(lzw, lzv, p*sizeof(*lzw)); CHKERRV(ierr);
	if(flowDesc->bc[0][ZPLUS].type != PERIODIC)
		lzw[p-1]--;
	if(numY == mesh->ny-1)
		lyw[n-1]++;
	numX = mesh->nx;
	numY = mesh->ny;
	numZ = (flowDesc->bc[0][ZPLUS].type == PERIODIC)? mesh->nz : mesh->nz-1;	
	ierr = DMDACreate3d(PETSC_COMM_WORLD, bx, by, bz, DMDA_STENCIL_BOX, numX, numY, numZ, m, n, p, 1, 1, lxw, lyw, lzw, &(wda)); CHKERRV(ierr);
	ierr = DMCompositeAddDM(pack, wda); CHKERRV(ierr);
	
	PetscFree(lxv);
	PetscFree(lyv);
	PetscFree(lzv);
	PetscFree(lxw);
	PetscFree(lyw);
	PetscFree(lzw);
		
	// create vectors
	// velocity fluxes
	ierr = DMCreateLocalVector(uda, &qxLocal); CHKERRV(ierr);
	ierr = DMCreateLocalVector(vda, &qyLocal); CHKERRV(ierr);
	ierr = DMCreateLocalVector(wda, &qzLocal); CHKERRV(ierr);

	// convection terms
	//ierr = DMCreateGlobalVector(uda, &(rxGlobal)); CHKERRV(ierr);	
	//ierr = DMCreateGlobalVector(vda, &(ryGlobal)); CHKERRV(ierr);
	//ierr = DMCreateGlobalVector(wda, &(rzGlobal)); CHKERRV(ierr);
}
