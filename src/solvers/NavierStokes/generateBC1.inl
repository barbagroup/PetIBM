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
	ierr = DMCompositeGetAccess(qPack, bc1, &bc1xGlobal, &bc1yGlobal); CHKERRV(ierr);
	               
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
					case DIRICHLET :
					case CONVECTIVE: bc1x[j][0] += coeffMinus*qx[j][-1]/mesh->dy[j]; break;
					default        : break;
				}
			}
			// +X
			if(mstart+m-1 == M-1) // if the +X face is in the current process
			{			
				switch(flowDesc->bc[0][XPLUS].type)
				{
					case DIRICHLET :
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
					case DIRICHLET :
					case CONVECTIVE: bc1x[0][i] += coeffMinus*qx[-1][i]; break;
					default        : break;
				}
			}
			// +Y
			if(nstart+n-1 == N-1) // if the +Y boundary is in the current process
			{
				switch(flowDesc->bc[0][YPLUS].type)
				{
					case DIRICHLET :
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
					case DIRICHLET :
					case CONVECTIVE: bc1y[j][0] += coeffMinus*qy[j][-1]; break;
					default        : break;
				}
			}
			// +X
			if(mstart+m-1 == M-1)
			{
				switch(flowDesc->bc[1][XPLUS].type)
				{
					case DIRICHLET :
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
					case DIRICHLET :
					case CONVECTIVE: bc1y[0][i] += coeffMinus*qy[-1][i]/mesh->dx[i]; break;
					default        : break;
				}
			}
			// +Y
			if(nstart+n-1 == N-1)
			{
				switch(flowDesc->bc[1][YPLUS].type)
				{
					case DIRICHLET :
					case CONVECTIVE: bc1y[N-1][i] += coeffPlus*qy[N][i]/mesh->dx[i]; break;
					default        : break;
				}
			}
		}
	}
	ierr = DMDAVecRestoreArray(vda, bc1yGlobal, &bc1y); CHKERRV(ierr);
	ierr = DMDAVecRestoreArray(vda, qyLocal, &qy); CHKERRV(ierr);

	ierr = DMCompositeRestoreAccess(qPack, bc1,  &bc1xGlobal, &bc1yGlobal); CHKERRV(ierr);
}

template <>
void NavierStokesSolver<3>::generateBC1()
{
	PetscErrorCode ierr;
	PetscInt       mstart, nstart, pstart, m, n, p, i, j, k, M, N, P;
	PetscReal      ***qx, ***qy, ***qz;
	PetscReal      ***bc1x, ***bc1y, ***bc1z;
	Vec            bc1xGlobal, bc1yGlobal, bc1zGlobal;
	PetscReal      nu = flowDesc->nu;
	PetscReal      alphaImplicit = simParams->alphaImplicit;
	PetscReal      coeffMinus = 0.0, coeffPlus = 0.0;

	ierr = VecSet(bc1, 0.0); CHKERRV(ierr);
	ierr = DMCompositeGetAccess(qPack, bc1, &bc1xGlobal, &bc1yGlobal, &bc1zGlobal); CHKERRV(ierr);
	
	// U-FLUXES
	ierr = DMDAVecGetArray(uda, bc1xGlobal, &bc1x); CHKERRV(ierr);
	ierr = DMDAVecGetArray(uda, qxLocal, &qx); CHKERRV(ierr);
	ierr = DMDAGetCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRV(ierr);
	ierr = DMDAGetInfo(uda, NULL, &M, &N, &P, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRV(ierr);
	// x-faces
	if(flowDesc->bc[0][XPLUS].type != PERIODIC) // don't update if the BC type is periodic
	{
		coeffMinus = alphaImplicit*nu*2.0/dxU[0]/(dxU[0]+dxU[1]);
		coeffPlus  = alphaImplicit*nu*2.0/dxU[M]/(dxU[M]+dxU[M-1]);
		//loop over all points on the x-face
		for(k=pstart; k<pstart+p; k++)
		{
			for(j=nstart; j<nstart+n; j++)
			{
				// -X
				if(mstart == 0) // if the -X face is in the current process
				{
					switch(flowDesc->bc[0][XMINUS].type)
					{
						case DIRICHLET :
						case CONVECTIVE: bc1x[k][j][0] += coeffMinus*qx[k][j][-1]/(mesh->dy[j]*mesh->dz[k]); break;
						default        : break;
					}
				}
				// +X
				if(mstart+m-1 == M-1) // if the +X face is in the current process
				{
					switch(flowDesc->bc[0][XPLUS].type)
					{
						case DIRICHLET :
						case CONVECTIVE: bc1x[k][j][M-1] += coeffPlus*qx[k][j][M]/(mesh->dy[j]*mesh->dz[k]); break;
						default        : break;
					}
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
		for(k=pstart; k<pstart+p; k++)
		{
			for(i=mstart; i<mstart+m; i++)
			{
				// -Y
				if(nstart == 0) // if the -Y boundary is in the current process
				{
					switch(flowDesc->bc[0][YMINUS].type)
					{
						case DIRICHLET :
						case CONVECTIVE: bc1x[k][0][i] += coeffMinus*qx[k][-1][i]; break;
						default        : break;
					}
				}
				// +Y
				if(nstart+n-1 == N-1) // if the +Y boundary is in the current process
				{
					switch(flowDesc->bc[0][YPLUS].type)
					{
						case DIRICHLET :
						case CONVECTIVE: bc1x[k][N-1][i] += coeffPlus*qx[k][N][i]; break;
						default        : break;
					}
				}
			}
		}
	}
	// z-faces
	if(flowDesc->bc[0][ZPLUS].type != PERIODIC) // don't update if the BC type is periodic
	{
		coeffMinus = alphaImplicit*nu*2.0/dzU[0]/(dzU[0]+dzU[1]);
		coeffPlus  = alphaImplicit*nu*2.0/dzU[P]/(dzU[P]+dzU[P-1]);
		// loop over all points on the z-face
		for(j=nstart; j<nstart+n; j++)
		{
			for(i=mstart; i<mstart+m; i++)
			{
				// -Z
				if(pstart == 0) // if the -Z boundary is in the current process
				{
					switch(flowDesc->bc[0][ZMINUS].type)
					{
						case DIRICHLET :
						case CONVECTIVE: bc1x[0][j][i] += coeffMinus*qx[-1][j][i]; break;
						default        : break;
					}
				}
				// +Z
				if(pstart+p-1 == P-1) // if the +Z boundary is in the current process
				{
					switch(flowDesc->bc[0][ZPLUS].type)
					{
						case DIRICHLET :
						case CONVECTIVE: bc1x[P-1][j][i] += coeffPlus*qx[P][j][i]; break;
						default        : break;
					}
				}
			}
		}
	}
	ierr = DMDAVecRestoreArray(uda, bc1xGlobal, &bc1x); CHKERRV(ierr);
	ierr = DMDAVecRestoreArray(uda, qxLocal, &qx); CHKERRV(ierr);
	
	// V-FLUXES
	ierr = DMDAVecGetArray(vda, bc1yGlobal, &bc1y); CHKERRV(ierr);
	ierr = DMDAVecGetArray(vda, qyLocal, &qy); CHKERRV(ierr);
	ierr = DMDAGetCorners(vda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRV(ierr);
	ierr = DMDAGetInfo(vda, NULL, &M, &N, &P, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRV(ierr);
	// x-faces
	if(flowDesc->bc[1][XPLUS].type != PERIODIC) // don't update if the BC type is periodic
	{
		coeffMinus = alphaImplicit*nu*2.0/dxV[0]/(dxV[0]+dxV[1]);
		coeffPlus  = alphaImplicit*nu*2.0/dxV[M]/(dxV[M]+dxV[M-1]);
		// loop over all points on the x-face
		for(k=pstart; k<pstart+p; k++)
		{
			for(j=nstart; j<nstart+n; j++)
			{
				// -X
				if(mstart == 0) // if the -X face is in the current process
				{
					switch(flowDesc->bc[1][XMINUS].type)
					{
						case DIRICHLET :
						case CONVECTIVE: bc1y[k][j][0] += coeffMinus*qy[k][j][-1]; break;
						default        : break;
					}
				}
				// +X
				if(mstart+m-1 == M-1) // if the +X face is in the current process
				{
					switch(flowDesc->bc[1][XPLUS].type)
					{
						case DIRICHLET :
						case CONVECTIVE: bc1y[k][j][M-1] += coeffPlus*qy[k][j][M]; break;
						default        : break;
					}
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
		for(k=pstart; k<pstart+p; k++)
		{
			for(i=mstart; i<mstart+m; i++)
			{
				// -Y
				if(nstart == 0) // if the -Y face is in the current process
				{
					switch(flowDesc->bc[1][YMINUS].type)
					{
						case DIRICHLET :
						case CONVECTIVE: bc1y[k][0][i] += coeffMinus*qy[k][-1][i]/(mesh->dz[k]*mesh->dx[i]); break;
						default        : break;
					}
				}
				// +Y
				if(nstart+n-1 == N-1) // if the +Y face is in the current process
				{
					switch(flowDesc->bc[1][YPLUS].type)
					{
						case DIRICHLET :
						case CONVECTIVE: bc1y[k][N-1][i] += coeffPlus*qy[k][N][i]/(mesh->dz[k]*mesh->dx[i]); break;
						default        : break;
					}
				}
			}
		}
	}
	// z-faces
	if(flowDesc->bc[1][ZPLUS].type != PERIODIC) // don't update if the BC type is periodic
	{
		coeffMinus = alphaImplicit*nu*2.0/dzV[0]/(dzV[0]+dzV[1]);
		coeffPlus  = alphaImplicit*nu*2.0/dzV[P]/(dzV[P]+dzV[P-1]);
		// loop over all points on the z-face
		for(j=nstart; j<nstart+n; j++)
		{
			for(i=mstart; i<mstart+m; i++)
			{
				// -Z
				if(pstart == 0) // if the -Z face is in the current process
				{
					switch(flowDesc->bc[1][ZMINUS].type)
					{
						case DIRICHLET :
						case CONVECTIVE: bc1y[0][j][i] += coeffMinus*qy[-1][j][i]; break;
						default        : break;
					}
				}
				// +Z
				if(pstart+p-1 == P-1) // if the +Z face is in the current process
				{
					switch(flowDesc->bc[1][ZPLUS].type)
					{
						case DIRICHLET :
						case CONVECTIVE: bc1y[P-1][j][i] += coeffPlus*qy[P][j][i]; break;
						default        : break;
					}
				}
			}
		}
	}
	ierr = DMDAVecRestoreArray(vda, bc1yGlobal, &bc1y); CHKERRV(ierr);
	ierr = DMDAVecRestoreArray(vda, qyLocal, &qy); CHKERRV(ierr);

	// W-FLUXES
	ierr = DMDAVecGetArray(wda, bc1zGlobal, &bc1z); CHKERRV(ierr);
	ierr = DMDAVecGetArray(wda, qzLocal, &qz); CHKERRV(ierr);
	ierr = DMDAGetCorners(wda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRV(ierr);
	ierr = DMDAGetInfo(wda, NULL, &M, &N, &P, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRV(ierr);
	// x-faces
	if(flowDesc->bc[2][XPLUS].type != PERIODIC) // don't update if the BC type is periodic
	{
		coeffMinus = alphaImplicit*nu*2.0/dxW[0]/(dxW[0]+dxW[1]);
		coeffPlus  = alphaImplicit*nu*2.0/dxW[M]/(dxW[M]+dxW[M-1]);
		// loop over all points on the x-face
		for(k=pstart; k<pstart+p; k++)
		{
			for(j=nstart; j<nstart+n; j++)
			{
				// -X
				if(mstart == 0) // if the -X face is in the current process
				{
					switch(flowDesc->bc[2][XMINUS].type)
					{
						case DIRICHLET :
						case CONVECTIVE: bc1z[k][j][0] += coeffMinus*qz[k][j][-1]; break;
						default        : break;
					}
				}
				// +X
				if(mstart+m-1 == M-1) // if the +X face is in the current process
				{
					switch(flowDesc->bc[2][XPLUS].type)
					{
						case DIRICHLET :
						case CONVECTIVE: bc1z[k][j][M-1] += coeffPlus*qz[k][j][M]; break;
						default        : break;
					}
				}
			}
		}
	}
	// y-faces
	if(flowDesc->bc[2][YPLUS].type != PERIODIC) // don't update if the BC type is periodic
	{
		coeffMinus = alphaImplicit*nu*2.0/dyW[0]/(dyW[0]+dyW[1]);
		coeffPlus  = alphaImplicit*nu*2.0/dyW[N]/(dyW[N]+dyW[N-1]);
		// loop over all points on the y-face
		for(k=pstart; k<pstart+p; k++)
		{
			for(i=mstart; i<mstart+m; i++)
			{
				// -Y
				if(nstart == 0) // if the -Y face is in the current process
				{
					switch(flowDesc->bc[2][YMINUS].type)
					{
						case DIRICHLET :
						case CONVECTIVE: bc1z[k][0][i] += coeffMinus*qz[k][-1][i]; break;
						default        : break;
					}
				}
				// +Y
				if(nstart+n-1 == N-1) // if the +Y face is in the current process
				{
					switch(flowDesc->bc[2][YPLUS].type)
					{
						case DIRICHLET :
						case CONVECTIVE: bc1z[k][N-1][i] += coeffPlus*qz[k][N][i]; break;
						default        : break;
					}
				}
			}
		}
	}
	// z-faces
	if(flowDesc->bc[2][ZPLUS].type != PERIODIC) // don't update if the BC type is periodic
	{
		coeffMinus = alphaImplicit*nu*2.0/dzW[0]/(dzW[0]+dzW[1]);
		coeffPlus  = alphaImplicit*nu*2.0/dzW[P]/(dzW[P]+dzW[P-1]);
		// loop over all points on the z-face
		for(j=nstart; j<nstart+n; j++)
		{
			for(i=mstart; i<mstart+m; i++)
			{
				// -Z
				if(pstart == 0) // if the -Z face is in the current process
				{
					switch(flowDesc->bc[2][ZMINUS].type)
					{
						case DIRICHLET :
						case CONVECTIVE: bc1z[0][j][i] += coeffMinus*qz[-1][j][i]/(mesh->dx[i]*mesh->dy[j]); break;
						default        : break;
					}
				}
				// +Z
				if(pstart+p-1 == P-1) // if the +Z face is in the current process
				{
					switch(flowDesc->bc[2][ZPLUS].type)
					{
						case DIRICHLET :
						case CONVECTIVE: bc1z[P-1][j][i] += coeffPlus*qz[P][j][i]/(mesh->dx[i]*mesh->dy[j]); break;
						default        : break;
					}
				}
			}
		}
	}
	ierr = DMDAVecRestoreArray(wda, bc1zGlobal, &bc1z); CHKERRV(ierr);
	ierr = DMDAVecRestoreArray(wda, qzLocal, &qz); CHKERRV(ierr);

	ierr = DMCompositeRestoreAccess(qPack, bc1, &bc1xGlobal, &bc1yGlobal, &bc1zGlobal); CHKERRV(ierr);
}
