/***************************************************************************//**
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
}

template <>
void NavierStokesSolver<2>::initializeMeshSpacings()
{
	PetscInt       numX, numY;
	PetscReal      dxMinus, dxPlus, dyMinus, dyPlus;

	// mesh spacings for U
	numX = (flowDesc->bc[0][XPLUS].type != PERIODIC)? mesh->nx-1 : mesh->nx; // number of U in the x-direction
	numY = mesh->ny; // number of U in the y-direction
	// dx
	dxU.resize(numX+1);
	for(PetscInt i=0; i<numX; i++)
	{
		dxU[i]   = mesh->dx[i];
		dxU[i+1] = (i < mesh->nx-1)? mesh->dx[i+1] : mesh->dx[0];
	}
	// dy
	dyU.resize(numY+1);
	for(PetscInt j=0; j<numY; j++)
	{
		// first check if the point is at the -Y or +Y edge of the mesh
		// then check if the boundary condition is periodic or not
		dyMinus = (j > 0)?          mesh->dy[j-1] : ((flowDesc->bc[0][YMINUS].type!=PERIODIC)? 0.0 : mesh->dy[mesh->ny-1]);
		dyPlus  = (j < mesh->ny-1)? mesh->dy[j+1] : ((flowDesc->bc[0][YPLUS].type !=PERIODIC)? 0.0 : mesh->dy[0]);

		dyU[j]   = 0.5*(mesh->dy[j] + dyMinus);
		dyU[j+1] = 0.5*(mesh->dy[j] + dyPlus);
	}

	// mesh spacings for V
	numX = mesh->nx; // number of V in the x-direction
	numY = (flowDesc->bc[1][YPLUS].type != PERIODIC)? mesh->ny-1 : mesh->ny; // number of V in the y-direction
	//dx
	dxV.resize(numX+1);
	for(PetscInt i=0; i<numX; i++)
	{
		// first check if the point is at the -X or +X edge of the mesh
		// then check if the boundary condition is periodic or not
		dxMinus = (i > 0)?          mesh->dx[i-1] : ((flowDesc->bc[1][XMINUS].type!=PERIODIC)? 0.0 : mesh->dx[mesh->nx-1]);
		dxPlus  = (i < mesh->nx-1)? mesh->dx[i+1] : ((flowDesc->bc[1][XPLUS].type!=PERIODIC)?  0.0 : mesh->dx[0]);

		dxV[i]   = 0.5*(mesh->dx[i] + dxMinus);
		dxV[i+1] = 0.5*(mesh->dx[i] + dxPlus);
	}
	//dy
	dyV.resize(numY+1);
	for(PetscInt j=0; j<numY; j++)
	{
		dyV[j]   = mesh->dy[j];
		dyV[j+1] = (j < mesh->ny-1)? mesh->dy[j+1] : mesh->dy[0];
	}
}

template <>
void NavierStokesSolver<3>::initializeMeshSpacings()
{
	PetscInt       numX, numY, numZ;
	PetscReal      dxMinus, dxPlus, dyMinus, dyPlus, dzMinus, dzPlus;

	// mesh spacings for U
	numX = (flowDesc->bc[0][XPLUS].type != PERIODIC)? mesh->nx-1 : mesh->nx; // number of U in the x-direction
	numY = mesh->ny; // number of U in the y-direction
	numZ = mesh->nz; // number of U in the z-direction
	// dx
	dxU.resize(numX+1);
	for(PetscInt i=0; i<numX; i++)
	{
		dxU[i]   = mesh->dx[i];
		dxU[i+1] = (i < mesh->nx-1)? mesh->dx[i+1] : mesh->dx[0];
	}
	// dy
	dyU.resize(numY+1);
	for(PetscInt j=0; j<numY; j++)
	{
		// first check if the point is at the -Y or +Y edge of the mesh
		// then check if the boundary condition is periodic or not
		dyMinus = (j > 0)?          mesh->dy[j-1] : ((flowDesc->bc[0][YMINUS].type!=PERIODIC)? 0.0 : mesh->dy[mesh->ny-1]);
		dyPlus  = (j < mesh->ny-1)? mesh->dy[j+1] : ((flowDesc->bc[0][YPLUS].type !=PERIODIC)? 0.0 : mesh->dy[0]);

		dyU[j]   = 0.5*(mesh->dy[j] + dyMinus);
		dyU[j+1] = 0.5*(mesh->dy[j] + dyPlus);
	}
	// dz
	dzU.resize(numZ+1);
	for(PetscInt k=0; k<numZ; k++)
	{
		// first check if the point is at the -Z or +Z edge of the mesh
		// then check if the boundary condition is periodic or not
		dzMinus = (k > 0)?          mesh->dz[k-1] : ((flowDesc->bc[0][ZMINUS].type!=PERIODIC)? 0.0 : mesh->dz[mesh->nz-1]);
		dzPlus  = (k < mesh->nz-1)? mesh->dz[k+1] : ((flowDesc->bc[0][ZPLUS].type !=PERIODIC)? 0.0 : mesh->dz[0]);

		dzU[k]   = 0.5*(mesh->dz[k] + dzMinus);
		dzU[k+1] = 0.5*(mesh->dz[k] + dzPlus);
	}

	// mesh spacings for V
	numX = mesh->nx; // number of V in the x-direction
	numY = (flowDesc->bc[1][YPLUS].type != PERIODIC)? mesh->ny-1 : mesh->ny; // number of V in the y-direction
	numZ = mesh->nz; // number of V in the z-direction
	//dx
	dxV.resize(numX+1);
	for(PetscInt i=0; i<numX; i++)
	{
		// first check if the point is at the -X or +X edge of the mesh
		// then check if the boundary condition is periodic or not
		dxMinus = (i > 0)?          mesh->dx[i-1] : ((flowDesc->bc[1][XMINUS].type!=PERIODIC)? 0.0 : mesh->dx[mesh->nx-1]);
		dxPlus  = (i < mesh->nx-1)? mesh->dx[i+1] : ((flowDesc->bc[1][XPLUS].type!=PERIODIC)?  0.0 : mesh->dx[0]);

		dxV[i]   = 0.5*(mesh->dx[i] + dxMinus);
		dxV[i+1] = 0.5*(mesh->dx[i] + dxPlus);
	}
	//dy
	dyV.resize(numY+1);
	for(PetscInt j=0; j<numY; j++)
	{
		dyV[j]   = mesh->dy[j];
		dyV[j+1] = (j < mesh->ny-1)? mesh->dy[j+1] : mesh->dy[0];
	}
	// dz
	dzV.resize(numZ+1);
	for(PetscInt k=0; k<numZ; k++)
	{
		// first check if the point is at the -Z or +Z edge of the mesh
		// then check if the boundary condition is periodic or not
		dzMinus = (k > 0)?          mesh->dz[k-1] : ((flowDesc->bc[1][ZMINUS].type!=PERIODIC)? 0.0 : mesh->dz[mesh->nz-1]);
		dzPlus  = (k < mesh->nz-1)? mesh->dz[k+1] : ((flowDesc->bc[1][ZPLUS].type !=PERIODIC)? 0.0 : mesh->dz[0]);

		dzV[k]   = 0.5*(mesh->dz[k] + dzMinus);
		dzV[k+1] = 0.5*(mesh->dz[k] + dzPlus);
	}

	// mesh spacings for W
	numX = mesh->nx; // number of W in the x-direction
	numY = mesh->ny; // number of W in the y-direction
	numZ = (flowDesc->bc[2][ZPLUS].type != PERIODIC)? mesh->nz-1 : mesh->nz; // number of W in the z-direction
	// dx
	dxW.resize(numX+1);
	for(PetscInt i=0; i<numX; i++)
	{
		// first check if the point is at the -X or +X edge of the mesh
		// then check if the boundary condition is periodic or not
		dxMinus = (i > 0)?          mesh->dx[i-1] : ((flowDesc->bc[2][XMINUS].type!=PERIODIC)? 0.0 : mesh->dx[mesh->nx-1]);
		dxPlus  = (i < mesh->nx-1)? mesh->dx[i+1] : ((flowDesc->bc[2][XPLUS].type!=PERIODIC)?  0.0 : mesh->dx[0]);

		dxW[i]   = 0.5*(mesh->dx[i] + dxMinus);
		dxW[i+1] = 0.5*(mesh->dx[i] + dxPlus);
	}
	// dy
	dyW.resize(numY+1);
	for(PetscInt j=0; j<numY; j++)
	{
		// first check if the point is at the -Y or +Y edge of the mesh
		// then check if the boundary condition is periodic or not
		dyMinus = (j > 0)?          mesh->dy[j-1] : ((flowDesc->bc[2][YMINUS].type!=PERIODIC)? 0.0 : mesh->dy[mesh->ny-1]);
		dyPlus  = (j < mesh->ny-1)? mesh->dy[j+1] : ((flowDesc->bc[2][YPLUS].type !=PERIODIC)? 0.0 : mesh->dy[0]);

		dyW[j]   = 0.5*(mesh->dy[j] + dyMinus);
		dyW[j+1] = 0.5*(mesh->dy[j] + dyPlus);
	}
	// dz
	dzW.resize(numZ+1);
	for(PetscInt k=0; k<numZ; k++)
	{
		dzW[k]   = mesh->dz[k];
		dzW[k+1] = (k < mesh->nz-1)? mesh->dz[k+1] : mesh->dz[0];
	}
}