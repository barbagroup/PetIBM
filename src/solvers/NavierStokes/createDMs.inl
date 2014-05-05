#include <petscdmcomposite.h>

template <>
void NavierStokesSolver<2>::createDMs()
{
	PetscErrorCode    ierr;
	PetscInt          m, n;
	const PetscInt    *lxp, *lyp;
	PetscInt          *lxu, *lyu, *lxv, *lyv;
	PetscInt          numX, numY;
	DMDABoundaryType  bx, by;
	
	// set boundary types for velocity flux variables
	bx = (flowDesc->bc[0][XPLUS].type == PERIODIC)? DMDA_BOUNDARY_PERIODIC : DMDA_BOUNDARY_GHOSTED;
	by = (flowDesc->bc[0][YPLUS].type == PERIODIC)? DMDA_BOUNDARY_PERIODIC : DMDA_BOUNDARY_GHOSTED;	
	// Create distributed array data structures
	// pressure
	numX = mesh->nx;
	numY = mesh->ny;
	ierr = DMDACreate2d(PETSC_COMM_WORLD, bx, by, DMDA_STENCIL_STAR, numX, numY, PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL, &pda); CHKERRV(ierr);
	ierr = DMDAGetOwnershipRanges(pda, &lxp, &lyp, NULL); CHKERRV(ierr);
	ierr = DMDAGetInfo(pda, NULL, NULL, NULL, NULL, &m, &n, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRV(ierr);
	// packed fluxes
	ierr = DMCompositeCreate(PETSC_COMM_WORLD, &pack); CHKERRV(ierr);
	// x-velocity
	ierr = PetscMalloc(m*sizeof(*lxu), &lxu); CHKERRV(ierr);
	ierr = PetscMalloc(n*sizeof(*lyu), &lyu); CHKERRV(ierr);
	ierr = PetscMemcpy(lxu, lxp, m*sizeof(*lxu)); CHKERRV(ierr);
	ierr = PetscMemcpy(lyu, lyp, n*sizeof(*lyu)); CHKERRV(ierr);
	numX = mesh->nx;
	numY = mesh->ny;
	if(flowDesc->bc[0][XPLUS].type != PERIODIC)
	{
		lxu[m-1]--;
		numX = mesh->nx-1;
	}
	ierr = DMDACreate2d(PETSC_COMM_WORLD, bx, by, DMDA_STENCIL_BOX, numX, numY, m, n, 1, 1, lxu, lyu, &uda); CHKERRV(ierr);
	ierr = DMCompositeAddDM(pack, uda); CHKERRV(ierr);
	// y-velocity
	ierr = PetscMalloc(m*sizeof(*lxv), &lxv); CHKERRV(ierr);
	ierr = PetscMalloc(n*sizeof(*lyv), &lyv); CHKERRV(ierr);
	ierr = PetscMemcpy(lxv, lxp, m*sizeof(*lxv)); CHKERRV(ierr);
	ierr = PetscMemcpy(lyv, lyp, n*sizeof(*lyv)); CHKERRV(ierr);
	numX = mesh->nx;
	numY = mesh->ny;
	if(flowDesc->bc[1][YPLUS].type != PERIODIC)
	{
		lyv[n-1]--;
		numY = mesh->ny-1;
	}
	ierr = DMDACreate2d(PETSC_COMM_WORLD, bx, by, DMDA_STENCIL_BOX, numX, numY, m, n, 1, 1, lxv, lyv, &vda); CHKERRV(ierr);
	ierr = DMCompositeAddDM(pack, vda); CHKERRV(ierr);
	
	PetscFree(lxu);
	PetscFree(lyu);	
	PetscFree(lxv);
	PetscFree(lyv);
}

template <>
void NavierStokesSolver<3>::createDMs()
{
	PetscErrorCode    ierr;
	PetscInt          m, n, p;
	const PetscInt    *lxp, *lyp, *lzp;
	PetscInt          *lxu, *lyu, *lzu;
	PetscInt          *lxv, *lyv, *lzv;
	PetscInt          *lxw, *lyw, *lzw;
	PetscInt          numX, numY, numZ;
	DMDABoundaryType  bx, by, bz;
	
	// set boundary types for velocity variables
	bx = (flowDesc->bc[0][XPLUS].type == PERIODIC)? DMDA_BOUNDARY_PERIODIC : DMDA_BOUNDARY_GHOSTED;
	by = (flowDesc->bc[0][YPLUS].type == PERIODIC)? DMDA_BOUNDARY_PERIODIC : DMDA_BOUNDARY_GHOSTED;
	bz = (flowDesc->bc[0][ZPLUS].type == PERIODIC)? DMDA_BOUNDARY_PERIODIC : DMDA_BOUNDARY_GHOSTED;
	// Create distributed array data structures
	// pressure
	numX = mesh->nx;
	numY = mesh->ny;
	numZ = mesh->nz;
	ierr = DMDACreate3d(PETSC_COMM_WORLD, bx, by, bz, DMDA_STENCIL_STAR, numX, numY, numZ, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL, NULL, &(uda)); CHKERRV(ierr);
	ierr = DMDAGetOwnershipRanges(pda, &lxp, &lyp, &lzp); CHKERRV(ierr);
	ierr = DMDAGetInfo(pda, NULL, NULL, NULL, NULL, &m, &n, &p, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRV(ierr);
	// packed fluxes
	ierr = DMCompositeCreate(PETSC_COMM_WORLD, &pack); CHKERRV(ierr);
	// x-velocity
	ierr = PetscMalloc(m*sizeof(*lxu), &lxu); CHKERRV(ierr);
	ierr = PetscMalloc(n*sizeof(*lyu), &lyu); CHKERRV(ierr);
	ierr = PetscMalloc(p*sizeof(*lzu), &lzu); CHKERRV(ierr);
	ierr = PetscMemcpy(lxu, lxp, m*sizeof(*lxu)); CHKERRV(ierr);
	ierr = PetscMemcpy(lyu, lyp, n*sizeof(*lyu)); CHKERRV(ierr);
	ierr = PetscMemcpy(lzu, lzp, p*sizeof(*lzu)); CHKERRV(ierr);
	numX = mesh->nx;
	numY = mesh->ny;
	numZ = mesh->nz;
	if(flowDesc->bc[0][XPLUS].type != PERIODIC)
	{
		lxu[m-1]--;
		numX = mesh->nx-1;
	}
	ierr = DMDACreate3d(PETSC_COMM_WORLD, bx, by, bz, DMDA_STENCIL_BOX, numX, numY, numZ, m, n, p, 1, 1, lxu, lyu, lzu, &uda); CHKERRV(ierr);	
	ierr = DMCompositeAddDM(pack, uda); CHKERRV(ierr);
	// y-velocity
	ierr = PetscMalloc(m*sizeof(*lxv), &lxv); CHKERRV(ierr);
	ierr = PetscMalloc(n*sizeof(*lyv), &lyv); CHKERRV(ierr);
	ierr = PetscMalloc(p*sizeof(*lzv), &lzv); CHKERRV(ierr);
	ierr = PetscMemcpy(lxv, lxp, m*sizeof(*lxv)); CHKERRV(ierr);
	ierr = PetscMemcpy(lyv, lyp, n*sizeof(*lyv)); CHKERRV(ierr);
	ierr = PetscMemcpy(lzv, lzp, p*sizeof(*lzv)); CHKERRV(ierr);
	numX = mesh->nx;
	numY = mesh->ny;
	numZ = mesh->nz;
	if(flowDesc->bc[1][YPLUS].type != PERIODIC)
	{
		lyv[n-1]--;
		numY = mesh->ny-1;
	}
	ierr = DMDACreate3d(PETSC_COMM_WORLD, bx, by, bz, DMDA_STENCIL_BOX, numX, numY, numZ, m, n, p, 1, 1, lxv, lyv, lzv, &vda); CHKERRV(ierr);
	ierr = DMCompositeAddDM(pack, vda); CHKERRV(ierr);
	// z-velocity
	ierr = PetscMalloc(m*sizeof(*lxw), &lxw); CHKERRV(ierr);
	ierr = PetscMalloc(n*sizeof(*lyw), &lyw); CHKERRV(ierr);
	ierr = PetscMalloc(p*sizeof(*lzw), &lzw); CHKERRV(ierr);
	ierr = PetscMemcpy(lxw, lxp, m*sizeof(*lxw)); CHKERRV(ierr);
	ierr = PetscMemcpy(lyw, lyp, n*sizeof(*lyw)); CHKERRV(ierr);
	ierr = PetscMemcpy(lzw, lzp, p*sizeof(*lzw)); CHKERRV(ierr);
	numX = mesh->nx;
	numY = mesh->ny;
	numZ = mesh->nz;
	if(flowDesc->bc[2][ZPLUS].type != PERIODIC)
	{
		lzw[p-1]--;
		numZ = mesh->nz-1;
	}
	ierr = DMDACreate3d(PETSC_COMM_WORLD, bx, by, bz, DMDA_STENCIL_BOX, numX, numY, numZ, m, n, p, 1, 1, lxw, lyw, lzw, &wda); CHKERRV(ierr);
	ierr = DMCompositeAddDM(pack, wda); CHKERRV(ierr);
	
	PetscFree(lxu);
	PetscFree(lyu);
	PetscFree(lzu);
	PetscFree(lxv);
	PetscFree(lyv);
	PetscFree(lzv);
	PetscFree(lxw);
	PetscFree(lyw);
	PetscFree(lzw);
}
