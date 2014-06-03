template <>
void NavierStokesSolver<2>::initialiseFluxes()
{
	PetscErrorCode ierr;
	PetscInt       mstart, nstart, m, n, i, j;
	PetscReal      **qx, **qy;
	               
	// U-FLUXES
	ierr = DMDAVecGetArray(uda, qxLocal, &qx); CHKERRV(ierr);
	ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRV(ierr);
	// Set interior values for u-fluxes
	for(j=nstart; j<nstart+n; j++)
	{
		for(i=mstart; i<mstart+m; i++)
		{
			qx[j][i] = flowDesc->initialVelocity[0] * mesh->dy[j];
		}
	}
	ierr = DMDAVecRestoreArray(uda, qxLocal, &qx); CHKERRV(ierr);
	// Update interior ghost cells for u-fluxes
	ierr = DMDALocalToLocalBegin(uda, qxLocal, INSERT_VALUES, qxLocal); CHKERRV(ierr);
	ierr = DMDALocalToLocalEnd(uda, qxLocal, INSERT_VALUES, qxLocal); CHKERRV(ierr);
	
	// V-FLUXES
	ierr = DMDAVecGetArray(vda, qyLocal, &qy); CHKERRV(ierr);
	ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRV(ierr);
	// Set interior values for v-fluxes
	for(j=nstart; j<nstart+n; j++)
	{
		for(i=mstart; i<mstart+m; i++)
		{
			qy[j][i] = flowDesc->initialVelocity[1] * mesh->dx[i];
		}
	}
	ierr = DMDAVecRestoreArray(vda, qyLocal, &qy); CHKERRV(ierr);
	// Update interior ghost cells for v-fluxes
	ierr = DMDALocalToLocalBegin(vda, qyLocal, INSERT_VALUES, qyLocal); CHKERRV(ierr);
	ierr = DMDALocalToLocalEnd(vda, qyLocal, INSERT_VALUES, qyLocal); CHKERRV(ierr);
}

template <>
void NavierStokesSolver<3>::initialiseFluxes()
{
	PetscErrorCode ierr;
	PetscInt       mstart, nstart, pstart, m, n, p, i, j, k;
	PetscReal      ***qx, ***qy, ***qz;
	               
	// U-FLUXES
	ierr = DMDAVecGetArray(uda, qxLocal, &qx); CHKERRV(ierr);
	ierr = DMDAGetCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRV(ierr);
	// Set interior values for u-fluxes
	for(k=pstart; k<pstart+p; k++)
	{
		for(j=nstart; j<nstart+n; j++)
		{
			for(i=mstart; i<mstart+m; i++)
			{
				qx[k][j][i] = flowDesc->initialVelocity[0] * (mesh->dy[j]*mesh->dz[k]);
			}
		}
	}
	ierr = DMDAVecRestoreArray(uda, qxLocal, &qx); CHKERRV(ierr);
	// Update interior ghost cells for u-fluxes
	ierr = DMDALocalToLocalBegin(uda, qxLocal, INSERT_VALUES, qxLocal); CHKERRV(ierr);
	ierr = DMDALocalToLocalEnd(uda, qxLocal, INSERT_VALUES, qxLocal); CHKERRV(ierr);
	
	// V-FLUXES
	ierr = DMDAVecGetArray(vda, qyLocal, &qy); CHKERRV(ierr);
	ierr = DMDAGetCorners(vda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRV(ierr);
	// Set interior values for v-fluxes
	for(k=pstart; k<pstart+p; k++)
	{
		for(j=nstart; j<nstart+n; j++)
		{
			for(i=mstart; i<mstart+m; i++)
			{
				qy[k][j][i] = flowDesc->initialVelocity[1] * (mesh->dx[i]*mesh->dz[k]);
			}
		}
	}
	ierr = DMDAVecRestoreArray(vda, qyLocal, &qy); CHKERRV(ierr);
	// Update interior ghost cells for v-fluxes
	ierr = DMDALocalToLocalBegin(vda, qyLocal, INSERT_VALUES, qyLocal); CHKERRV(ierr);
	ierr = DMDALocalToLocalEnd(vda, qyLocal, INSERT_VALUES, qyLocal); CHKERRV(ierr);

	// W-FLUXES
	ierr = DMDAVecGetArray(wda, qzLocal, &qz); CHKERRV(ierr);
	ierr = DMDAGetCorners(wda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRV(ierr);
	// Set interior values for w-fluxes
	for(k=pstart; k<pstart+p; k++)
	{
		for(j=nstart; j<nstart+n; j++)
		{
			for(i=mstart; i<mstart+m; i++)
			{
				qz[k][j][i] = flowDesc->initialVelocity[2] * (mesh->dx[i]*mesh->dy[j]);
			}
		}
	}
	ierr = DMDAVecRestoreArray(wda, qzLocal, &qz); CHKERRV(ierr);
	// Update interior ghost cells for w-fluxes
	ierr = DMDALocalToLocalBegin(wda, qzLocal, INSERT_VALUES, qzLocal); CHKERRV(ierr);
	ierr = DMDALocalToLocalEnd(wda, qzLocal, INSERT_VALUES, qzLocal); CHKERRV(ierr);
}
