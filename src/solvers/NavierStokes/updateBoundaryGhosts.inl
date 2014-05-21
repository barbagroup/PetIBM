template <>
void NavierStokesSolver<2>::updateBoundaryGhosts()
{
	PetscErrorCode ierr;
	PetscInt       mstart, nstart, m, n, i, j, M, N;
	PetscReal      **qx, **qy;
	               
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
					case NEUMANN   : qx[j][-1] = qx[j][0]; break;
					default        : break;
				}
			}
			// +X
			if(mstart+m-1 == M-1) // if the +X face is in the current process
			{			
				switch(flowDesc->bc[0][XPLUS].type)
				{
					case DIRICHLET : qx[j][M] = flowDesc->bc[0][XPLUS].value*mesh->dy[j]; break;
					case NEUMANN   : qx[j][M] = qx[j][M-1]; break;
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
					case NEUMANN   : qx[-1][i] = qx[0][i]/mesh->dy[0]; break;
					default        : break;
				}
			}
			// +Y
			if(nstart+n-1 == N-1) // if the +Y boundary is in the current process
			{
				switch(flowDesc->bc[0][YPLUS].type)
				{
					case DIRICHLET : qx[N][i] = flowDesc->bc[0][YPLUS].value; break;
					case NEUMANN   : qx[N][i] = qx[N-1][i]/mesh->dy[N-1]; break;
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
					case NEUMANN   : qy[j][-1] = qy[j][0]/mesh->dx[0]; break;
					default        : break;
				}
			}
			// +X
			if(mstart+m-1 == M-1)
			{
				switch(flowDesc->bc[1][XPLUS].type)
				{
					case DIRICHLET : qy[j][M] = flowDesc->bc[1][XPLUS].value; break;
					case NEUMANN   : qy[j][M] = qy[j][M-1]/mesh->dx[M-1]; break;
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
					case NEUMANN   : qy[-1][i] = qy[0][i]; break;
					default        : break;
				}
			}
			// +Y
			if(nstart+n-1 == N-1)
			{
				switch(flowDesc->bc[1][YPLUS].type)
				{
					case DIRICHLET : qy[N][i] = flowDesc->bc[1][YPLUS].value*mesh->dx[i]; break;
					case NEUMANN   : qy[N][i] = qy[N-1][i]; break;
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
}
