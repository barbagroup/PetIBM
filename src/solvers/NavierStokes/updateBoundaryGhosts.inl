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
		for(j=nstart; j<nstart+n; j++)
		{
			// -X
			if(mstart == 0) // if the -X face is in the current process
			{
				switch(flowDesc->bc[0][XMINUS].type)
				{
					case DIRICHLET : qx[j][mstart-1] = flowDesc->bc[0][XMINUS].value*mesh->dy[j]; break;
					case NEUMANN   : qx[j][mstart-1] = qx[j][mstart]; break;
					default        : break;
				}
			}
			// +X
			if(mstart+m-1 == M-1) // if the +X face is in the current process
			{			
				switch(flowDesc->bc[0][XPLUS].type)
				{
					case DIRICHLET : qx[j][mstart+m] = flowDesc->bc[0][XPLUS].value*mesh->dy[j]; break;
					case NEUMANN   : qx[j][mstart+m] = qx[j][mstart+m-1]; break;
					default        : break;
				}
			}
		}
	}
	// y-faces
	if(flowDesc->bc[0][YPLUS].type != PERIODIC) // don't update if the BC type is periodic
	{
		for(i=mstart; i<mstart+m; i++)
		{	
			// -Y
			if(nstart == 0) // if the -Y boundary is in the current process
			{
				switch(flowDesc->bc[0][YMINUS].type)
				{
					case DIRICHLET : qx[nstart-1][i] = flowDesc->bc[0][YMINUS].value; break;
					case NEUMANN   : qx[nstart-1][i] = qx[nstart][i]/mesh->dy[nstart]; break;
					default        : break;
				}
			}
			// +Y
			if(nstart+n-1 == N-1) // if the +Y boundary is in the current process
			{
				switch(flowDesc->bc[0][YPLUS].type)
				{
					case DIRICHLET : qx[nstart+n][i] = flowDesc->bc[0][YPLUS].value; break;
					case NEUMANN   : qx[nstart+n][i] = qx[nstart+n-1][i]/mesh->dy[nstart+n-1]; break;
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
		for(j=nstart; j<nstart+n; j++)
		{
			// the following values are set for ghost cells in all processes
			// a LocalToLocal copy is called at the end of this function
			// to ensure that all the internal ghost values are correct
			
			// -X
			if(mstart == 0)
			{
				switch(flowDesc->bc[1][XMINUS].type)
				{
					case DIRICHLET : qy[j][mstart-1] = flowDesc->bc[1][XMINUS].value; break;
					case NEUMANN   : qy[j][mstart-1] = qy[j][mstart] / mesh->dx[mstart]; break;
					default        : break;
				}
			}
			// +X
			if(mstart+m-1 == M-1)
			{
				switch(flowDesc->bc[1][XPLUS].type)
				{
					case DIRICHLET : qy[j][mstart+m] = flowDesc->bc[1][XPLUS].value; break;
					case NEUMANN   : qy[j][mstart+m] = qy[j][mstart+m-1] / mesh->dx[mstart+m-1]; break;
					default        : break;
				}
			}
		}
	}
	// y-faces
	if(flowDesc->bc[1][YPLUS].type != PERIODIC) // don't update if the BC type is periodic
	{
		for(i=mstart; i<mstart+m; i++)
		{	
			// -Y
			if(nstart == 0)
			{
				switch(flowDesc->bc[1][YMINUS].type)
				{
					case DIRICHLET : qy[nstart-1][i] = flowDesc->bc[1][YMINUS].value*mesh->dx[i]; break;
					case NEUMANN   : qy[nstart-1][i] = qy[nstart][i]; break;
					default        : break;
				}
			}
			// +Y
			if(nstart+n-1 == N-1)
			{
				switch(flowDesc->bc[1][YPLUS].type)
				{
					case DIRICHLET : qy[nstart+n][i] = flowDesc->bc[1][YPLUS].value*mesh->dx[i]; break;
					case NEUMANN   : qy[nstart+n][i] = qy[nstart+n-1][i]; break;
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
