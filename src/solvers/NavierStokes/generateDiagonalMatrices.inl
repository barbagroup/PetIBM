template <>
PetscErrorCode NavierStokesSolver<2>::generateDiagonalMatrices()
{
	PetscErrorCode ierr;
	PetscInt       mstart, nstart, m, n, i, j;
	Vec            MHatxGlobal, MHatyGlobal;
	Vec            RInvxGlobal, RInvyGlobal;
	Vec            BNxGlobal, BNyGlobal;
	PetscReal      **MHatx, **MHaty;
	PetscReal      **RInvx, **RInvy;
	PetscReal      **BNx, **BNy;

	ierr = DMCompositeGetAccess(qPack, MHat, &MHatxGlobal, &MHatyGlobal); CHKERRQ(ierr);
	ierr = DMCompositeGetAccess(qPack, RInv, &RInvxGlobal, &RInvyGlobal); CHKERRQ(ierr);
	ierr = DMCompositeGetAccess(qPack, BN,   &BNxGlobal,   &BNyGlobal); CHKERRQ(ierr);

	// x-direction
	ierr = DMDAVecGetArray(uda, MHatxGlobal, &MHatx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(uda, RInvxGlobal, &RInvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(uda, BNxGlobal,   &BNx); CHKERRQ(ierr);
	ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
	for(j=nstart; j<nstart+n; j++)
	{
		for(i=mstart; i<mstart+m; i++)
		{
			MHatx[j][i] = (i<mesh->nx-1)? 0.5*(mesh->dx[i] + mesh->dx[i+1]) : 0.5*(mesh->dx[i] + mesh->dx[0]);
			RInvx[j][i] = 1.0/mesh->dy[j];
			BNx[j][i]   = simParams->dt/(MHatx[j][i]*RInvx[j][i]);
		}
	}
	ierr = DMDAVecRestoreArray(uda, MHatxGlobal, &MHatx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(uda, RInvxGlobal, &RInvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(uda, BNxGlobal,   &BNx); CHKERRQ(ierr);

	// y-direction
	ierr = DMDAVecGetArray(vda, MHatyGlobal, &MHaty); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(vda, RInvyGlobal, &RInvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(vda, BNyGlobal,   &BNy); CHKERRQ(ierr);
	ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
	for(j=nstart; j<nstart+n; j++)
	{
		for(i=mstart; i<mstart+m; i++)
		{
			MHaty[j][i] = (j<mesh->ny-1)? 0.5*(mesh->dy[j] + mesh->dy[j+1]) : 0.5*(mesh->dy[j] + mesh->dy[0]);
			RInvy[j][i] = 1.0/mesh->dx[i];
			BNy[j][i]   = simParams->dt/(MHaty[j][i]*RInvy[j][i]);
		}
	}
	ierr = DMDAVecRestoreArray(vda, MHatyGlobal, &MHaty); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(vda, RInvyGlobal, &RInvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(vda, BNyGlobal,   &BNy); CHKERRQ(ierr);

	ierr = DMCompositeRestoreAccess(qPack, MHat, &MHatxGlobal, &MHatyGlobal); CHKERRQ(ierr);
	ierr = DMCompositeRestoreAccess(qPack, RInv, &RInvxGlobal, &RInvyGlobal); CHKERRQ(ierr);
	ierr = DMCompositeRestoreAccess(qPack, BN,   &BNxGlobal,   &BNyGlobal); CHKERRQ(ierr);

	return 0;
}

template <>
PetscErrorCode NavierStokesSolver<3>::generateDiagonalMatrices()
{
	PetscErrorCode ierr;
	PetscInt       mstart, nstart, pstart, m, n, p, i, j, k;
	Vec            MHatxGlobal, MHatyGlobal, MHatzGlobal;
	Vec            RInvxGlobal, RInvyGlobal, RInvzGlobal;
	Vec            BNxGlobal, BNyGlobal, BNzGlobal;
	PetscReal      ***MHatx, ***MHaty, ***MHatz;
	PetscReal      ***RInvx, ***RInvy, ***RInvz;
	PetscReal      ***BNx, ***BNy, ***BNz;
	
	ierr = DMCompositeGetAccess(qPack, MHat, &MHatxGlobal, &MHatyGlobal, &MHatzGlobal); CHKERRQ(ierr);
	ierr = DMCompositeGetAccess(qPack, RInv, &RInvxGlobal, &RInvyGlobal, &RInvzGlobal); CHKERRQ(ierr);
	ierr = DMCompositeGetAccess(qPack, BN,   &BNxGlobal,   &BNyGlobal,   &BNzGlobal); CHKERRQ(ierr);
	
	// x-direction
	ierr = DMDAVecGetArray(uda, MHatxGlobal, &MHatx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(uda, RInvxGlobal, &RInvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(uda, BNxGlobal,   &BNx); CHKERRQ(ierr);
	ierr = DMDAGetCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
	for(k=pstart; k<pstart+p; k++)
	{
		for(j=nstart; j<nstart+n; j++)
		{
			for(i=mstart; i<mstart+m; i++)
			{
				MHatx[k][j][i] = (i < mesh->nx-1)? 0.5*(mesh->dx[i] + mesh->dx[i+1]) : 0.5*(mesh->dx[i] + mesh->dx[0]);
				RInvx[k][j][i] = 1.0/(mesh->dy[j]*mesh->dz[k]);
				BNx[k][j][i]   = simParams->dt/(MHatx[k][j][i]*RInvx[k][j][i]);
			}
		}
	}
	ierr = DMDAVecRestoreArray(uda, MHatxGlobal, &MHatx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(uda, RInvxGlobal, &RInvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(uda, BNxGlobal,   &BNx); CHKERRQ(ierr);
	
	// y-direction
	ierr = DMDAVecGetArray(vda, MHatyGlobal, &MHaty); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(vda, RInvyGlobal, &RInvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(vda, BNyGlobal,   &BNy); CHKERRQ(ierr);
	ierr = DMDAGetCorners(vda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
	for(k=pstart; k<pstart+p; k++)
	{
		for(j=nstart; j<nstart+n; j++)
		{
			for(i=mstart; i<mstart+m; i++)
			{
				MHaty[k][j][i] = (j < mesh->ny-1)? 0.5*(mesh->dy[j] + mesh->dy[j+1]) : 0.5*(mesh->dy[j] + mesh->dy[0]);
				RInvy[k][j][i] = 1.0/(mesh->dz[k]*mesh->dx[i]);
				BNy[k][j][i]   = simParams->dt/(MHaty[k][j][i]*RInvy[k][j][i]);
			}
		}
	}
	ierr = DMDAVecRestoreArray(vda, MHatyGlobal, &MHaty); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(vda, RInvyGlobal, &RInvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(vda, BNyGlobal,   &BNy); CHKERRQ(ierr);
	
	// z-direction
	ierr = DMDAVecGetArray(wda, MHatzGlobal, &MHatz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(wda, RInvzGlobal, &RInvz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(wda, BNzGlobal,   &BNz); CHKERRQ(ierr);
	ierr = DMDAGetCorners(wda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
	for(k=pstart; k<pstart+p; k++)
	{
		for(j=nstart; j<nstart+n; j++)
		{
			for(i=mstart; i<mstart+m; i++)
			{
				MHatz[k][j][i] = (k < mesh->nz-1)? 0.5*(mesh->dz[k] + mesh->dz[k+1]) : 0.5*(mesh->dz[k] + mesh->dz[0]);
				RInvz[k][j][i] = 1.0/(mesh->dx[i]*mesh->dy[j]);
				BNz[k][j][i]   = simParams->dt/(MHatz[k][j][i]*RInvz[k][j][i]);
			}
		}
	}
	ierr = DMDAVecRestoreArray(wda, MHatzGlobal, &MHatz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(wda, RInvzGlobal, &RInvz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(wda, BNzGlobal,   &BNz); CHKERRQ(ierr);
	
	ierr = DMCompositeRestoreAccess(qPack, MHat, &MHatxGlobal, &MHatyGlobal, &MHatzGlobal); CHKERRQ(ierr);
	ierr = DMCompositeRestoreAccess(qPack, RInv, &RInvxGlobal, &RInvyGlobal, &RInvzGlobal); CHKERRQ(ierr);
	ierr = DMCompositeRestoreAccess(qPack, BN,   &BNxGlobal,   &BNyGlobal,   &BNzGlobal); CHKERRQ(ierr);

	return 0;
}