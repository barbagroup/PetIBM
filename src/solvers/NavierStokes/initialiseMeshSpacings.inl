template<>
void NavierStokesSolver<2>::initialiseMeshSpacings()
{
	PetscInt         numX, numY;
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

template<>
void NavierStokesSolver<3>::initialiseMeshSpacings()
{
	
}