template <>
void NavierStokesSolver<2>::initialiseFluxes()
{
	PetscErrorCode ierr;
	PetscInt       mstart, nstart, m, n;
	PetscReal      **qx, **qy;
	Vec            qxGlobal, qyGlobal;
	
	ierr = DMCompositeGetAccess(qPack, q, &qxGlobal, &qyGlobal); CHKERRV(ierr);
	               
	// U-FLUXES
	ierr = DMDAVecGetArray(uda, qxGlobal, &qx); CHKERRV(ierr);
	ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRV(ierr);
	// Set interior values for u-fluxes
	for(PetscInt j=nstart; j<nstart+n; j++)
	{
		for(PetscInt i=mstart; i<mstart+m; i++)
		{
			qx[j][i] = flowDesc->initialVelocity[0] * mesh->dy[j];
		}
	}
	ierr = DMDAVecRestoreArray(uda, qxGlobal, &qx); CHKERRV(ierr);
	
	// V-FLUXES
	ierr = DMDAVecGetArray(vda, qyGlobal, &qy); CHKERRV(ierr);
	ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRV(ierr);
	// Set interior values for v-fluxes
	for(PetscInt j=nstart; j<nstart+n; j++)
	{
		for(PetscInt i=mstart; i<mstart+m; i++)
		{
			qy[j][i] = flowDesc->initialVelocity[1] * mesh->dx[i];
		}
	}
	ierr = DMDAVecRestoreArray(vda, qyGlobal, &qy); CHKERRV(ierr);
	
	ierr = DMCompositeRestoreAccess(qPack, q, &qxGlobal, &qyGlobal); CHKERRV(ierr);
	ierr = DMCompositeScatter(qPack, q, qxLocal, qyLocal); CHKERRV(ierr);
}

template <>
void NavierStokesSolver<3>::initialiseFluxes()
{
	PetscErrorCode ierr;
	PetscInt       mstart, nstart, pstart, m, n, p;
	PetscReal      ***qx, ***qy, ***qz;
	Vec            qxGlobal, qyGlobal, qzGlobal;
	
	ierr = DMCompositeGetAccess(qPack, q, &qxGlobal, &qyGlobal, &qzGlobal); CHKERRV(ierr);
	               
	// U-FLUXES
	ierr = DMDAVecGetArray(uda, qxGlobal, &qx); CHKERRV(ierr);
	ierr = DMDAGetCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRV(ierr);
	// Set interior values for u-fluxes
	for(PetscInt k=pstart; k<pstart+p; k++)
	{
		for(PetscInt j=nstart; j<nstart+n; j++)
		{
			for(PetscInt i=mstart; i<mstart+m; i++)
			{
				qx[k][j][i] = flowDesc->initialVelocity[0] * (mesh->dy[j]*mesh->dz[k]);
			}
		}
	}
	ierr = DMDAVecRestoreArray(uda, qxGlobal, &qx); CHKERRV(ierr);
	
	// V-FLUXES
	ierr = DMDAVecGetArray(vda, qyGlobal, &qy); CHKERRV(ierr);
	ierr = DMDAGetCorners(vda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRV(ierr);
	// Set interior values for v-fluxes
	for(PetscInt k=pstart; k<pstart+p; k++)
	{
		for(PetscInt j=nstart; j<nstart+n; j++)
		{
			for(PetscInt i=mstart; i<mstart+m; i++)
			{
				qy[k][j][i] = flowDesc->initialVelocity[1] * (mesh->dx[i]*mesh->dz[k]);
			}
		}
	}
	ierr = DMDAVecRestoreArray(vda, qyGlobal, &qy); CHKERRV(ierr);

	// W-FLUXES
	ierr = DMDAVecGetArray(wda, qzGlobal, &qz); CHKERRV(ierr);
	ierr = DMDAGetCorners(wda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRV(ierr);
	// Set interior values for w-fluxes
	for(PetscInt k=pstart; k<pstart+p; k++)
	{
		for(PetscInt j=nstart; j<nstart+n; j++)
		{
			for(PetscInt i=mstart; i<mstart+m; i++)
			{
				qz[k][j][i] = flowDesc->initialVelocity[2] * (mesh->dx[i]*mesh->dy[j]);
			}
		}
	}
	ierr = DMDAVecRestoreArray(wda, qzGlobal, &qz); CHKERRV(ierr);
	
	ierr = DMCompositeRestoreAccess(qPack, q, &qxGlobal, &qyGlobal, &qzGlobal); CHKERRV(ierr);
	ierr = DMCompositeScatter(qPack, q, qxLocal, qyLocal, qzLocal); CHKERRV(ierr);
}
