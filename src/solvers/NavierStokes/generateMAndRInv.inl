template <>
void NavierStokesSolver<2>::generateMAndRInv()
{
	PetscErrorCode ierr;
	PetscInt       mstart, nstart, m, n, i, j;
	Vec            MxGlobal, MyGlobal;
	Vec            RInvxGlobal, RInvyGlobal;
	PetscReal      **Mx, **My;
	PetscReal      **RInvx, **RInvy;

	ierr = DMCompositeGetAccess(pack, M, &MxGlobal, &MyGlobal); CHKERRV(ierr);
	ierr = DMCompositeGetAccess(pack, RInv, &RInvxGlobal, &RInvyGlobal); CHKERRV(ierr);

	// x-direction
	ierr = DMDAVecGetArray(uda, MxGlobal, &Mx); CHKERRV(ierr);
	ierr = DMDAVecGetArray(uda, RInvxGlobal, &RInvx); CHKERRV(ierr);
	ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRV(ierr);
	for(j=nstart; j<nstart+n; j++)
	{
		for(i=mstart; i<mstart+m; i++)
		{
			Mx[j][i] = (i<mesh->nx-1)? 0.5*(mesh->dx[i] + mesh->dx[i+1]) : 0.5*(mesh->dx[i] + mesh->dx[0]);
			RInvx[j][i] = 1.0/mesh->dy[j];
		}
	}
	ierr = DMDAVecRestoreArray(uda, MxGlobal, &Mx); CHKERRV(ierr);
	ierr = DMDAVecRestoreArray(uda, RInvxGlobal, &RInvx); CHKERRV(ierr);

	// y-direction
	ierr = DMDAVecGetArray(vda, MyGlobal, &My); CHKERRV(ierr);
	ierr = DMDAVecGetArray(vda, RInvyGlobal, &RInvy); CHKERRV(ierr);
	ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRV(ierr);
	for(j=nstart; j<nstart+n; j++)
	{
		for(i=mstart; i<mstart+m; i++)
		{
			My[j][i] = (j<mesh->ny-1)? 0.5*(mesh->dy[j] + mesh->dy[j+1]) : 0.5*(mesh->dy[j] + mesh->dy[0]);
			RInvy[j][i] = 1.0/mesh->dx[i];
		}
	}
	ierr = DMDAVecRestoreArray(vda, MyGlobal, &My); CHKERRV(ierr);
	ierr = DMDAVecRestoreArray(vda, RInvyGlobal, &RInvy); CHKERRV(ierr);

	ierr = DMCompositeRestoreAccess(pack, M, &MxGlobal, &MyGlobal); CHKERRV(ierr);
	ierr = DMCompositeRestoreAccess(pack, RInv, &RInvxGlobal, &RInvyGlobal); CHKERRV(ierr);
}

template <>
void NavierStokesSolver<3>::generateMAndRInv()
{
	PetscErrorCode ierr;
	PetscInt       mstart, nstart, pstart, m, n, p, i, j, k;
	Vec            MxGlobal, MyGlobal, MzGlobal;
	Vec            RInvxGlobal, RInvyGlobal, RInvzGlobal;
	PetscReal      ***Mx, ***My, ***Mz;
	PetscReal      ***RInvx, ***RInvy, ***RInvz;
	
	ierr = DMCompositeGetAccess(pack, M, &MxGlobal, &MyGlobal, &MzGlobal); CHKERRV(ierr);
	ierr = DMCompositeGetAccess(pack, RInv, &RInvxGlobal, &RInvyGlobal, &RInvzGlobal); CHKERRV(ierr);
	
	// x-direction
	ierr = DMDAVecGetArray(uda, MxGlobal, &Mx); CHKERRV(ierr);
	ierr = DMDAVecGetArray(uda, RInvxGlobal, &RInvx); CHKERRV(ierr);
	ierr = DMDAGetCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRV(ierr);
	for(k=pstart; k<pstart+p; k++)
	{
		for(j=nstart; j<nstart+n; j++)
		{
			for(i=mstart; i<mstart+m; i++)
			{
				Mx[k][j][i] = (i < mesh->nx-1)? 0.5*(mesh->dx[i] + mesh->dx[i+1]) : 0.5*(mesh->dx[i] + mesh->dx[0]);
				RInvx[k][j][i] = 1.0/(mesh->dy[j]*mesh->dz[k]);
			}
		}
	}
	ierr = DMDAVecRestoreArray(uda, MxGlobal, &Mx); CHKERRV(ierr);
	ierr = DMDAVecRestoreArray(uda, RInvxGlobal, &RInvx); CHKERRV(ierr);
	
	// y-direction
	ierr = DMDAVecGetArray(vda, MyGlobal, &My); CHKERRV(ierr);
	ierr = DMDAVecGetArray(vda, RInvyGlobal, &RInvy); CHKERRV(ierr);
	ierr = DMDAGetCorners(vda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRV(ierr);
	for(k=pstart; k<pstart+p; k++)
	{
		for(j=nstart; j<nstart+n; j++)
		{
			for(i=mstart; i<mstart+m; i++)
			{
				My[k][j][i] = (j < mesh->ny-1)? 0.5*(mesh->dy[j] + mesh->dy[j+1]) : 0.5*(mesh->dy[j] + mesh->dy[0]);
				RInvy[k][j][i] = 1.0/(mesh->dx[i]*mesh->dz[k]);
			}
		}
	}
	ierr = DMDAVecRestoreArray(vda, MyGlobal, &My); CHKERRV(ierr);
	ierr = DMDAVecRestoreArray(vda, RInvyGlobal, &RInvy); CHKERRV(ierr);
	
	// z-direction
	ierr = DMDAVecGetArray(wda, MzGlobal, &Mz); CHKERRV(ierr);
	ierr = DMDAVecGetArray(wda, RInvzGlobal, &RInvz); CHKERRV(ierr);
	ierr = DMDAGetCorners(wda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRV(ierr);
	for(k=pstart; k<pstart+p; k++)
	{
		for(j=nstart; j<nstart+n; j++)
		{
			for(i=mstart; i<mstart+m; i++)
			{
				Mz[k][j][i] = (k < mesh->nz-1)? 0.5*(mesh->dz[k] + mesh->dz[k+1]) : 0.5*(mesh->dz[k] + mesh->dz[0]);
				RInvz[k][j][i] = 1.0/(mesh->dy[j]*mesh->dz[k]);
			}
		}
	}
	ierr = DMDAVecRestoreArray(wda, MzGlobal, &Mz); CHKERRV(ierr);
	ierr = DMDAVecRestoreArray(wda, RInvzGlobal, &RInvz); CHKERRV(ierr);
	
	ierr = DMCompositeRestoreAccess(pack, M, &MxGlobal, &MyGlobal, &MzGlobal); CHKERRV(ierr);
	ierr = DMCompositeRestoreAccess(pack, RInv, &RInvxGlobal, &RInvyGlobal, &RInvzGlobal); CHKERRV(ierr);
}