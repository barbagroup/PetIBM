template <>
void NavierStokesSolver<2>::updateBoundaryGhosts()
{
	PetscErrorCode ierr;
	PetscInt       mstart, nstart, m, n, i, j, M, N;
	PetscReal      **qx, **qy;
	PetscReal      dt = simParams->dt;
	PetscReal      beta;
	               
	// U-FLUXES
	ierr = DMDAVecGetArray(uda, qxLocal, &qx); CHKERRV(ierr);
	ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRV(ierr);
	ierr = DMDAGetInfo(uda, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRV(ierr);
	// x-faces
	if(flowDesc->bc[0][XPLUS].type != PERIODIC) // don't update if the BC type is periodic
	{
		// loop over all points on the x-face
		for(j=nstart; j<nstart+n; j++)
		{
			// -X
			if(mstart == 0) // if the -X face is in the current process
			{
				switch(flowDesc->bc[0][XMINUS].type)
				{
					case DIRICHLET : qx[j][-1] = flowDesc->bc[0][XMINUS].value*mesh->dy[j]; break;
					case CONVECTIVE: beta = flowDesc->bc[0][XMINUS].value*dt/dxU[0];
					                 qx[j][-1] = (1-beta)*qx[j][-1] + beta*qx[j][0]; break;
					default        : break;
				}
			}
			// +X
			if(mstart+m-1 == M-1) // if the +X face is in the current process
			{			
				switch(flowDesc->bc[0][XPLUS].type)
				{
					case DIRICHLET : qx[j][M] = flowDesc->bc[0][XPLUS].value*mesh->dy[j]; break;
					case CONVECTIVE: beta = flowDesc->bc[0][XPLUS].value*dt/dxU[M];
					                 qx[j][M] = (1-beta)*qx[j][M] + beta*qx[j][M-1]; break;
					default        : break;
				}
			}
		}
	}
	// y-faces
	if(flowDesc->bc[0][YPLUS].type != PERIODIC) // don't update if the BC type is periodic
	{
		// loop over all points on the y-face
		for(i=mstart; i<mstart+m; i++)
		{	
			// -Y
			if(nstart == 0) // if the -Y boundary is in the current process
			{
				switch(flowDesc->bc[0][YMINUS].type)
				{
					case DIRICHLET : qx[-1][i] = flowDesc->bc[0][YMINUS].value; break;
					case CONVECTIVE: beta = flowDesc->bc[0][YMINUS].value*dt/dyU[0];
					                 qx[-1][i] = (1.0-beta)*qx[-1][i] + beta*qx[0][i]/mesh->dy[0]; break;
					default        : break;
				}
			}
			// +Y
			if(nstart+n-1 == N-1) // if the +Y boundary is in the current process
			{
				switch(flowDesc->bc[0][YPLUS].type)
				{
					case DIRICHLET : qx[N][i] = flowDesc->bc[0][YPLUS].value; break;
					case CONVECTIVE: beta = flowDesc->bc[0][YPLUS].value*dt/dyU[N];
					                 qx[N][i] = (1.0-beta)*qx[N][i] + beta*qx[N-1][i]/mesh->dy[N-1]; break;
					default        : break;
				}
			}
		}
	}
	ierr = DMDAVecRestoreArray(uda, qxLocal, &qx); CHKERRV(ierr);
		
	
	// V-FLUXES
	ierr = DMDAVecGetArray(vda, qyLocal, &qy); CHKERRV(ierr);
	ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRV(ierr);
	ierr = DMDAGetInfo(vda, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRV(ierr);
	// x-faces
	if(flowDesc->bc[1][XPLUS].type != PERIODIC) // don't update if the BC type is periodic
	{
		// loop over all points on the x-face
		for(j=nstart; j<nstart+n; j++)
		{
			// -X
			if(mstart == 0)
			{
				switch(flowDesc->bc[1][XMINUS].type)
				{
					case DIRICHLET : qy[j][-1] = flowDesc->bc[1][XMINUS].value; break;
					case CONVECTIVE: beta = flowDesc->bc[1][XMINUS].value*dt/dxV[0];
					                 qy[j][-1] = (1.0-beta)*qy[j][-1] + beta*qy[j][0]/mesh->dx[0]; break;
					default        : break;
				}
			}
			// +X
			if(mstart+m-1 == M-1)
			{
				switch(flowDesc->bc[1][XPLUS].type)
				{
					case DIRICHLET : qy[j][M] = flowDesc->bc[1][XPLUS].value; break;
					case CONVECTIVE: beta = flowDesc->bc[1][XPLUS].value*dt/dxV[M];
					                 qy[j][M] = (1.0-beta)*qy[j][M] + beta*qy[j][M-1]/mesh->dx[M-1]; break;
					default        : break;
				}
			}
		}
	}
	// y-faces
	if(flowDesc->bc[1][YPLUS].type != PERIODIC) // don't update if the BC type is periodic
	{
		// loop over all points on the y-face
		for(i=mstart; i<mstart+m; i++)
		{	
			// -Y
			if(nstart == 0)
			{
				switch(flowDesc->bc[1][YMINUS].type)
				{
					case DIRICHLET : qy[-1][i] = flowDesc->bc[1][YMINUS].value*mesh->dx[i]; break;
					case CONVECTIVE: beta = flowDesc->bc[1][YMINUS].value*dt/dyV[0];
					                 qy[-1][i] = (1.0-beta)*qy[-1][i] + beta*qy[0][i]; break;
					default        : break;
				}
			}
			// +Y
			if(nstart+n-1 == N-1)
			{
				switch(flowDesc->bc[1][YPLUS].type)
				{
					case DIRICHLET : qy[N][i] = flowDesc->bc[1][YPLUS].value*mesh->dx[i]; break;
					case CONVECTIVE: beta = flowDesc->bc[1][YPLUS].value*dt/dyV[N];
					                 qy[N][i] = (1.0-beta)*qy[N][i] + beta*qy[N-1][i]; break;
					default        : break;
				}
			}
		}
	}
	ierr = DMDAVecRestoreArray(vda, qyLocal, &qy); CHKERRV(ierr);
}

template <>
void NavierStokesSolver<3>::updateBoundaryGhosts()
{
	PetscErrorCode ierr;
	PetscInt       mstart, nstart, pstart, m, n, p, i, j, k, M, N, P;
	PetscReal      ***qx, ***qy, ***qz;
	PetscReal      dt = simParams->dt;
	PetscReal      beta;

	// U-FLUXES
	ierr = DMDAVecGetArray(uda, qxLocal, &qx); CHKERRV(ierr);
	ierr = DMDAGetCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRV(ierr);
	ierr = DMDAGetInfo(uda, NULL, &M, &N, &P, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRV(ierr);
	// x-faces
	if(flowDesc->bc[0][XPLUS].type != PERIODIC) // don't update if the BC type is periodic
	{
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
						case DIRICHLET : qx[k][j][-1] = flowDesc->bc[0][XMINUS].value*(mesh->dy[j]*mesh->dz[k]); break;
						case CONVECTIVE: beta = flowDesc->bc[0][XMINUS].value*dt/dxU[0];
					                     qx[k][j][-1] = (1-beta)*qx[k][j][-1] + beta*qx[k][j][0]; break;
					    default        : break;
					}
				}
				// +X
				if(mstart+m-1 == M-1) // if the +X face is in the current process
				{
					switch(flowDesc->bc[0][XPLUS].type)
					{
						case DIRICHLET : qx[k][j][M] = flowDesc->bc[0][XPLUS].value*(mesh->dy[j]*mesh->dz[k]); break;
						case CONVECTIVE: beta = flowDesc->bc[0][XPLUS].value*dt/dxU[M];
					                     qx[k][j][M] = (1-beta)*qx[k][j][M] + beta*qx[k][j][M-1]; break;
						default        : break;
					}
				}
			}
		}
	}
	// y-faces
	if(flowDesc->bc[0][YPLUS].type != PERIODIC) // don't update if the BC type is periodic
	{
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
						case DIRICHLET : qx[k][-1][i] = flowDesc->bc[0][YMINUS].value; break;
						case CONVECTIVE: beta = flowDesc->bc[0][YMINUS].value*dt/dyU[0];
						                 qx[k][-1][i] = (1.0-beta)*qx[k][-1][i] + beta*qx[k][0][i]/(mesh->dy[0]*mesh->dz[k]); break;
						default        : break;
					}
				}
				// +Y
				if(nstart+n-1 == N-1) // if the +Y boundary is in the current process
				{
					switch(flowDesc->bc[0][YPLUS].type)
					{
						case DIRICHLET : qx[k][N][i] = flowDesc->bc[0][YPLUS].value; break;
						case CONVECTIVE: beta = flowDesc->bc[0][YPLUS].value*dt/dyU[N];
						                 qx[k][N][i] = (1.0-beta)*qx[k][N][i] + beta*qx[k][N-1][i]/(mesh->dy[N-1]*mesh->dz[k]); break;
						default        : break;
					}
				}
			}
		}
	}
	// z-faces
	if(flowDesc->bc[0][ZPLUS].type != PERIODIC) // don't update if the BC type is periodic
	{
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
						case DIRICHLET : qx[-1][j][i] = flowDesc->bc[0][ZMINUS].value; break;
						case CONVECTIVE: beta = flowDesc->bc[0][ZMINUS].value*dt/dzU[0];
						                 qx[-1][j][i] = (1.0-beta)*qx[-1][j][i] + beta*qx[0][j][i]/(mesh->dy[j]*mesh->dz[0]); break;
						default        : break;
					}
				}
				// +Z
				if(pstart+p-1 == P-1) // if the +Z boundary is in the current process
				{
					switch(flowDesc->bc[0][ZPLUS].type)
					{
						case DIRICHLET : qx[P][j][i] = flowDesc->bc[0][ZPLUS].value; break;
						case CONVECTIVE: beta = flowDesc->bc[0][ZPLUS].value*dt/dzU[P];
						                 qx[P][j][i] = (1.0-beta)*qx[P][j][i] + beta*qx[P-1][j][i]/(mesh->dy[j]*mesh->dz[P-1]); break;
						default        : break;
					}
				}
			}
		}
	}
	ierr = DMDAVecRestoreArray(uda, qxLocal, &qx); CHKERRV(ierr);
	
	// V-FLUXES
	ierr = DMDAVecGetArray(vda, qyLocal, &qy); CHKERRV(ierr);
	ierr = DMDAGetCorners(vda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRV(ierr);
	ierr = DMDAGetInfo(vda, NULL, &M, &N, &P, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRV(ierr);
	// x-faces
	if(flowDesc->bc[1][XPLUS].type != PERIODIC) // don't update if the BC type is periodic
	{
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
						case DIRICHLET : qy[k][j][-1] = flowDesc->bc[1][XMINUS].value; break;
						case CONVECTIVE: beta = flowDesc->bc[1][XMINUS].value*dt/dxV[0];
						                 qy[k][j][-1] = (1.0-beta)*qy[k][j][-1] + beta*qy[k][j][0]/(mesh->dx[0]*mesh->dz[k]); break;
						default        : break;
					}
				}
				// +X
				if(mstart+m-1 == M-1) // if the +X face is in the current process
				{
					switch(flowDesc->bc[1][XPLUS].type)
					{
						case DIRICHLET : qy[k][j][M] = flowDesc->bc[1][XPLUS].value; break;
						case CONVECTIVE: beta = flowDesc->bc[1][XPLUS].value*dt/dxV[M];
						                 qy[k][j][M] = (1.0-beta)*qy[k][j][M] + beta*qy[k][j][M-1]/(mesh->dx[M-1]*mesh->dz[k]); break;
						default        : break;
					}
				}
			}
		}
	}
	// y-faces
	if(flowDesc->bc[1][YPLUS].type != PERIODIC) // don't update if the BC type is periodic
	{
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
						case DIRICHLET : qy[k][-1][i] = flowDesc->bc[1][YMINUS].value*(mesh->dz[k]*mesh->dx[i]); break;
						case CONVECTIVE: beta = flowDesc->bc[1][YMINUS].value*dt/dyV[0];
						                 qy[k][-1][i] = (1.0-beta)*qy[k][-1][i] + beta*qy[k][0][i]; break;
						default        : break;
					}
				}
				// +Y
				if(nstart+n-1 == N-1) // if the +Y face is in the current process
				{
					switch(flowDesc->bc[1][YPLUS].type)
					{
						case DIRICHLET : qy[k][N][i] = flowDesc->bc[1][YPLUS].value*(mesh->dz[k]*mesh->dx[i]); break;
						case CONVECTIVE: beta = flowDesc->bc[1][YPLUS].value*dt/dyV[N];
						                 qy[k][N][i] = (1.0-beta)*qy[k][N][i] + beta*qy[k][N-1][i]; break;
						default        : break;
					}
				}
			}
		}
	}
	// z-faces
	if(flowDesc->bc[1][ZPLUS].type != PERIODIC) // don't update if the BC type is periodic
	{
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
						case DIRICHLET : qy[-1][j][i] = flowDesc->bc[1][ZMINUS].value; break;
						case CONVECTIVE: beta = flowDesc->bc[1][ZMINUS].value*dt/dzV[0];
						                 qy[-1][j][i] = (1.0-beta)*qy[-1][j][i] + beta*qy[0][j][i]/(mesh->dx[i]*mesh->dz[0]); break;
						default        : break;
					}
				}
				// +Z
				if(pstart+p-1 == P-1) // if the +Z face is in the current process
				{
					switch(flowDesc->bc[1][ZPLUS].type)
					{
						case DIRICHLET : qy[P][j][i] = flowDesc->bc[1][ZPLUS].value; break;
						case CONVECTIVE: beta = flowDesc->bc[1][ZPLUS].value*dt/dzV[P];
						                 qy[P][j][i] = (1.0-beta)*qy[P][j][i] + beta*qy[P-1][j][i]/(mesh->dx[i]*mesh->dz[P-1]); break;
						default        : break;
					}
				}
			}
		}
	}
	ierr = DMDAVecRestoreArray(vda, qyLocal, &qy); CHKERRV(ierr);

	// W-FLUXES
	ierr = DMDAVecGetArray(wda, qzLocal, &qz); CHKERRV(ierr);
	ierr = DMDAGetCorners(wda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRV(ierr);
	ierr = DMDAGetInfo(wda, NULL, &M, &N, &P, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRV(ierr);
	// x-faces
	if(flowDesc->bc[2][XPLUS].type != PERIODIC) // don't update if the BC type is periodic
	{
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
						case DIRICHLET : qz[k][j][-1] = flowDesc->bc[2][XMINUS].value; break;
						case CONVECTIVE: beta = flowDesc->bc[2][XMINUS].value*dt/dxW[0];
						                 qz[k][j][-1] = (1.0-beta)*qz[k][j][-1] + beta*qz[k][j][0]/(mesh->dx[0]*mesh->dy[j]); break;
						default        : break;
					}
				}
				// +X
				if(mstart+m-1 == M-1) // if the +X face is in the current process
				{
					switch(flowDesc->bc[2][XPLUS].type)
					{
						case DIRICHLET : qz[k][j][M] = flowDesc->bc[2][XPLUS].value; break;
						case CONVECTIVE: beta = flowDesc->bc[2][XPLUS].value*dt/dxW[M];
						                 qz[k][j][M] = (1.0-beta)*qz[k][j][M] + beta*qz[k][j][M-1]/(mesh->dx[M-1]*mesh->dy[j]); break;
						default        : break;
					}
				}
			}
		}
	}
	// y-faces
	if(flowDesc->bc[2][YPLUS].type != PERIODIC) // don't update if the BC type is periodic
	{
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
						case DIRICHLET : qz[k][-1][i] = flowDesc->bc[2][YMINUS].value; break;
						case CONVECTIVE: beta = flowDesc->bc[2][YMINUS].value*dt/dyW[0];
						                 qz[k][-1][i] = (1.0-beta)*qz[k][-1][i] + beta*qz[k][0][i]/(mesh->dx[i]*mesh->dy[0]); break;
						default        : break;
					}
				}
				// +Y
				if(nstart+n-1 == N-1) // if the +Y face is in the current process
				{
					switch(flowDesc->bc[2][YPLUS].type)
					{
						case DIRICHLET : qz[k][N][i] = flowDesc->bc[2][YPLUS].value; break;
						case CONVECTIVE: beta = flowDesc->bc[2][YPLUS].value*dt/dyW[N];
						                 qz[k][N][i] = (1.0-beta)*qz[k][N][i] + beta*qz[k][N-1][i]/(mesh->dx[i]*mesh->dy[N-1]); break;
						default        : break;
					}
				}
			}
		}
	}
	// z-faces
	if(flowDesc->bc[2][ZPLUS].type != PERIODIC) // don't update if the BC type is periodic
	{
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
						case DIRICHLET : qz[-1][j][i] = flowDesc->bc[2][ZMINUS].value*(mesh->dx[i]*mesh->dy[j]); break;
						case CONVECTIVE: beta = flowDesc->bc[2][ZMINUS].value*dt/dzW[0];
						                 qz[-1][j][i] = (1.0-beta)*qz[-1][j][i] + beta*qz[0][j][i]; break;
						default        : break;
					}
				}
				// +Z
				if(pstart+p-1 == P-1) // if the +Z face is in the current process
				{
					switch(flowDesc->bc[2][ZPLUS].type)
					{
						case DIRICHLET : qz[P][j][i] = flowDesc->bc[2][ZPLUS].value*(mesh->dx[i]*mesh->dy[j]); break;
						case CONVECTIVE: beta = flowDesc->bc[2][ZPLUS].value*dt/dzW[P];
						                 qz[P][j][i] = (1.0-beta)*qz[P][j][i] + beta*qz[P-1][j][i]; break;
						default        : break;
					}
				}
			}
		}
	}
	ierr = DMDAVecRestoreArray(wda, qzLocal, &qz); CHKERRV(ierr);
}
