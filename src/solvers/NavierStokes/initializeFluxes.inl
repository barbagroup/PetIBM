template <>
PetscErrorCode NavierStokesSolver<2>::initializeFluxes()
{
	PetscErrorCode ierr;
	Vec            qxGlobal, qyGlobal;
	
	ierr = DMCompositeGetAccess(qPack, q, &qxGlobal, &qyGlobal); CHKERRQ(ierr);
	
	if(simParams->restart)
	{
		ierr = readFluxes(qxGlobal, qyGlobal); CHKERRQ(ierr);
	}
	else
	{
		PetscInt  mstart, nstart, m, n;
		PetscReal **qx, **qy;
		PetscReal initVel[2]  = {flowDesc->initialVelocity[0], flowDesc->initialVelocity[1]};
		PetscReal initPert[2] = {flowDesc->initialPerturbation[0], flowDesc->initialPerturbation[1]};
		PetscReal width[2]    = {mesh->x[mesh->nx] - mesh->x[0], mesh->y[mesh->ny] - mesh->y[0]};

		// U-FLUXES
		ierr = DMDAVecGetArray(uda, qxGlobal, &qx); CHKERRQ(ierr);
		ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
		// Set interior values for u-fluxes
		for(PetscInt j=nstart; j<nstart+n; j++)
		{
			for(PetscInt i=mstart; i<mstart+m; i++)
			{
				PetscReal x = -PETSC_PI + 2*PETSC_PI*(mesh->x[i+1] - mesh->x[0])/width[0],
				          y = -PETSC_PI + 2*PETSC_PI*(0.5*(mesh->y[j]+mesh->y[j+1]) - mesh->y[0])/width[1];
				
				qx[j][i] = (initVel[0] + initPert[0]*cos(0.5*x)*sin(y))*mesh->dy[j];
			}
		}
		ierr = DMDAVecRestoreArray(uda, qxGlobal, &qx); CHKERRQ(ierr);
		
		// V-FLUXES
		ierr = DMDAVecGetArray(vda, qyGlobal, &qy); CHKERRQ(ierr);
		ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
		// Set interior values for v-fluxes
		for(PetscInt j=nstart; j<nstart+n; j++)
		{
			for(PetscInt i=mstart; i<mstart+m; i++)
			{
				PetscReal x = -PETSC_PI + 2*PETSC_PI*(0.5*(mesh->x[i]+mesh->x[i+1]) - mesh->x[0])/width[0],
				          y = -PETSC_PI + 2*PETSC_PI*(mesh->y[j+1] - mesh->y[0])/width[1];
				
				qy[j][i] = (initVel[1] + initPert[1]*sin(x)*cos(0.5*y))*mesh->dx[i];
			}
		}
		ierr = DMDAVecRestoreArray(vda, qyGlobal, &qy); CHKERRQ(ierr);
	}
	
	ierr = DMCompositeRestoreAccess(qPack, q, &qxGlobal, &qyGlobal); CHKERRQ(ierr);
	
	ierr = DMCompositeScatter(qPack, q, qxLocal, qyLocal); CHKERRQ(ierr);

	return 0;
}

template <>
PetscErrorCode NavierStokesSolver<3>::initializeFluxes()
{
	PetscErrorCode ierr;
	Vec            qxGlobal, qyGlobal, qzGlobal;
	
	ierr = DMCompositeGetAccess(qPack, q, &qxGlobal, &qyGlobal, &qzGlobal); CHKERRQ(ierr);

	if(simParams->restart)
	{
		ierr = readFluxes(qxGlobal, qyGlobal, qzGlobal); CHKERRQ(ierr);
	}
	else
	{
		PetscInt  mstart, nstart, pstart, m, n, p;
		PetscReal ***qx, ***qy, ***qz;
		PetscReal initVel[3]  = {flowDesc->initialVelocity[0], flowDesc->initialVelocity[1], flowDesc->initialVelocity[2]};
		PetscReal initPert[3] = {flowDesc->initialPerturbation[0], flowDesc->initialPerturbation[1], flowDesc->initialPerturbation[2]};
		PetscReal width[3]    = {mesh->x[mesh->nx] - mesh->x[0], mesh->y[mesh->ny] - mesh->y[0], mesh->z[mesh->nz] - mesh->z[0]};
		
		// U-FLUXES
		ierr = DMDAVecGetArray(uda, qxGlobal, &qx); CHKERRQ(ierr);
		ierr = DMDAGetCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
		// Set interior values for u-fluxes
		for(PetscInt k=pstart; k<pstart+p; k++)
		{
			for(PetscInt j=nstart; j<nstart+n; j++)
			{
				for(PetscInt i=mstart; i<mstart+m; i++)
				{
					PetscReal x = -PETSC_PI + 2*PETSC_PI*(mesh->x[i+1] - mesh->x[0])/width[0],
				              y = -PETSC_PI + 2*PETSC_PI*(0.5*(mesh->y[j]+mesh->y[j+1]) - mesh->y[0])/width[1],
					          z = -PETSC_PI + 2*PETSC_PI*(0.5*(mesh->z[k]+mesh->z[k+1]) - mesh->z[0])/width[2];
					
					qx[k][j][i] = (initVel[0] + initPert[0]*cos(0.5*x)*sin(y)*sin(z))*(mesh->dy[j]*mesh->dz[k]);
				}
			}
		}
		ierr = DMDAVecRestoreArray(uda, qxGlobal, &qx); CHKERRQ(ierr);
		
		// V-FLUXES
		ierr = DMDAVecGetArray(vda, qyGlobal, &qy); CHKERRQ(ierr);
		ierr = DMDAGetCorners(vda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
		// Set interior values for v-fluxes
		for(PetscInt k=pstart; k<pstart+p; k++)
		{
			for(PetscInt j=nstart; j<nstart+n; j++)
			{
				for(PetscInt i=mstart; i<mstart+m; i++)
				{
					PetscReal x = -PETSC_PI + 2*PETSC_PI*(0.5*(mesh->x[i]+mesh->x[i+1]) - mesh->x[0])/width[0],
				              y = -PETSC_PI + 2*PETSC_PI*(mesh->y[j+1] - mesh->y[0])/width[1],
				              z = -PETSC_PI + 2*PETSC_PI*(0.5*(mesh->z[k]+mesh->z[k+1]) - mesh->z[0])/width[2];
				
					qy[k][j][i] = (initVel[1] + initPert[1]*sin(x)*cos(0.5*y)*sin(z))*(mesh->dx[i]*mesh->dz[k]);
				}
			}
		}
		ierr = DMDAVecRestoreArray(vda, qyGlobal, &qy); CHKERRQ(ierr);

		// W-FLUXES
		ierr = DMDAVecGetArray(wda, qzGlobal, &qz); CHKERRQ(ierr);
		ierr = DMDAGetCorners(wda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
		// Set interior values for w-fluxes
		for(PetscInt k=pstart; k<pstart+p; k++)
		{
			for(PetscInt j=nstart; j<nstart+n; j++)
			{
				for(PetscInt i=mstart; i<mstart+m; i++)
				{
					PetscReal x = -PETSC_PI + 2*PETSC_PI*(0.5*(mesh->x[i]+mesh->x[i+1]) - mesh->x[0])/width[0],
					          y = -PETSC_PI + 2*PETSC_PI*(0.5*(mesh->y[j]+mesh->y[j+1]) - mesh->y[0])/width[1],
					          z = -PETSC_PI + 2*PETSC_PI*(mesh->z[k+1] - mesh->z[0])/width[2];
					
					qz[k][j][i] = (initVel[2] + initPert[2]*sin(x)*sin(y)*cos(0.5*z))*(mesh->dx[i]*mesh->dy[j]);
				}
			}
		}
		ierr = DMDAVecRestoreArray(wda, qzGlobal, &qz); CHKERRQ(ierr);
	}
	
	ierr = DMCompositeRestoreAccess(qPack, q, &qxGlobal, &qyGlobal, &qzGlobal); CHKERRQ(ierr);
	
	ierr = DMCompositeScatter(qPack, q, qxLocal, qyLocal, qzLocal); CHKERRQ(ierr);

	return 0;
}
