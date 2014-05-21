template <>
void NavierStokesSolver<2>::generateBC1()
{
	PetscErrorCode ierr;
	PetscInt       mstart, nstart, m, n, i, j, M, N;
	PetscReal      **qx, **qy;
	PetscReal      **bc1x, **bc1y;
	Vec            bc1xGlobal, bc1yGlobal;
	PetscReal      nu = flowDesc->nu;
	PetscReal      alphaImplicit = simParams->alphaImplicit;
	PetscReal      coeffMinus = 0.0, coeffPlus = 0.0;

	ierr = VecSet(bc1, 0.0); CHKERRV(ierr);
	ierr = DMCompositeGetAccess(pack, bc1,  &bc1xGlobal, &bc1yGlobal); CHKERRV(ierr);
	               
	// U-FLUXES
	ierr = DMDAVecGetArray(uda, bc1xGlobal, &bc1x); CHKERRV(ierr);
	ierr = DMDAVecGetArray(uda, qxLocal, &qx); CHKERRV(ierr);
	ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRV(ierr);
	ierr = DMDAGetInfo(uda, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRV(ierr);
	// x-faces
	if(flowDesc->bc[0][XPLUS].type != PERIODIC) // don't update if the BC type is periodic
	{
		coeffMinus = alphaImplicit*nu*2.0/dxU[0]/(dxU[0]+dxU[1]);
		coeffPlus  = alphaImplicit*nu*2.0/dxU[M]/(dxU[M]+dxU[M-1]);
		//loop over all points on the x-face
		for(j=nstart; j<nstart+n; j++)
		{
			// -X
			if(mstart == 0) // if the -X face is in the current process
			{
				switch(flowDesc->bc[0][XMINUS].type)
				{
					case DIRICHLET : bc1x[j][0] += coeffMinus*qx[j][-1]/mesh->dy[j]; break;
					case CONVECTIVE: bc1x[j][0] += coeffMinus*qx[j][-1]/mesh->dy[j]; break;
					default        : break;
				}
			}
			// +X
			if(mstart+m-1 == M-1) // if the +X face is in the current process
			{			
				switch(flowDesc->bc[0][XPLUS].type)
				{
					case DIRICHLET : bc1x[j][M-1] += coeffPlus*qx[j][M]/mesh->dy[j]; break;
					case CONVECTIVE: bc1x[j][M-1] += coeffPlus*qx[j][M]/mesh->dy[j]; break;
					default        : break;
				}
			}
		}
	}
	// y-faces
	if(flowDesc->bc[0][YPLUS].type != PERIODIC) // don't update if the BC type is periodic
	{
		coeffMinus = alphaImplicit*nu*2.0/dyU[0]/(dyU[0]+dyU[1]);
		coeffPlus  = alphaImplicit*nu*2.0/dyU[N]/(dyU[N]+dyU[N-1]);
		// loop over all points on the y-face
		for(i=mstart; i<mstart+m; i++)
		{	
			// -Y
			if(nstart == 0) // if the -Y boundary is in the current process
			{
				switch(flowDesc->bc[0][YMINUS].type)
				{
					case DIRICHLET : bc1x[0][i] += coeffMinus*qx[-1][i]; break;
					case CONVECTIVE: bc1x[0][i] += coeffMinus*qx[-1][i]; break;
					default        : break;
				}
			}
			// +Y
			if(nstart+n-1 == N-1) // if the +Y boundary is in the current process
			{
				switch(flowDesc->bc[0][YPLUS].type)
				{
					case DIRICHLET : bc1x[N-1][i] += coeffPlus*qx[N][i]; break;
					case CONVECTIVE: bc1x[N-1][i] += coeffPlus*qx[N][i]; break;
					default        : break;
				}
			}
		}
	}
	ierr = DMDAVecRestoreArray(uda, bc1xGlobal, &bc1x); CHKERRV(ierr);
	ierr = DMDAVecRestoreArray(uda, qxLocal, &qx); CHKERRV(ierr);
		
	
	// V-FLUXES
	ierr = DMDAVecGetArray(vda, bc1yGlobal, &bc1y); CHKERRV(ierr);
	ierr = DMDAVecGetArray(vda, qyLocal, &qy); CHKERRV(ierr);
	ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRV(ierr);
	ierr = DMDAGetInfo(vda, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRV(ierr);
	// x-faces
	if(flowDesc->bc[1][XPLUS].type != PERIODIC) // don't update if the BC type is periodic
	{
		coeffMinus = alphaImplicit*nu*2.0/dxV[0]/(dxV[0]+dxV[1]);
		coeffPlus  = alphaImplicit*nu*2.0/dxV[M]/(dxV[M]+dxV[M-1]);
		// loop over all points on the x-face
		for(j=nstart; j<nstart+n; j++)
		{
			// -X
			if(mstart == 0)
			{
				switch(flowDesc->bc[1][XMINUS].type)
				{
					case DIRICHLET : bc1y[j][0] += coeffMinus*qy[j][-1]; break;
					case CONVECTIVE: bc1y[j][0] += coeffMinus*qy[j][-1]; break;
					default        : break;
				}
			}
			// +X
			if(mstart+m-1 == M-1)
			{
				switch(flowDesc->bc[1][XPLUS].type)
				{
					case DIRICHLET : bc1y[j][M-1] += coeffPlus*qy[j][M]; break;
					case CONVECTIVE: bc1y[j][M-1] += coeffPlus*qy[j][M]; break;
					default        : break;
				}
			}
		}
	}
	// y-faces
	if(flowDesc->bc[1][YPLUS].type != PERIODIC) // don't update if the BC type is periodic
	{
		coeffMinus = alphaImplicit*nu*2.0/dyV[0]/(dyV[0]+dyV[1]);
		coeffPlus  = alphaImplicit*nu*2.0/dyV[N]/(dyV[N]+dyV[N-1]);
		// loop over all points on the y-face
		for(i=mstart; i<mstart+m; i++)
		{	
			// -Y
			if(nstart == 0)
			{
				switch(flowDesc->bc[1][YMINUS].type)
				{
					case DIRICHLET : bc1y[0][i] += coeffMinus*qy[-1][i]/mesh->dx[i]; break;
					case CONVECTIVE: bc1y[0][i] += coeffMinus*qy[-1][i]/mesh->dx[i]; break;
					default        : break;
				}
			}
			// +Y
			if(nstart+n-1 == N-1)
			{
				switch(flowDesc->bc[1][YPLUS].type)
				{
					case DIRICHLET : bc1y[N-1][i] += coeffPlus*qy[N][i]/mesh->dx[i]; break;
					case CONVECTIVE: bc1y[N-1][i] += coeffPlus*qy[N][i]/mesh->dx[i]; break;
					default        : break;
				}
			}
		}
	}
	ierr = DMDAVecRestoreArray(vda, bc1yGlobal, &bc1y); CHKERRV(ierr);
	ierr = DMDAVecRestoreArray(vda, qyLocal, &qy); CHKERRV(ierr);

	ierr = DMCompositeRestoreAccess(pack, bc1,  &bc1xGlobal, &bc1yGlobal); CHKERRV(ierr);
}

template <>
void NavierStokesSolver<3>::generateBC1()
{
}
