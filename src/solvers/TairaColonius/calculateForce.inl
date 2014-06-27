template <>
PetscErrorCode TairaColoniusSolver<2>::calculateForce()
{
	PetscErrorCode ierr;
	PetscInt       mstart, nstart, m, n;
	Vec            fGlobal, fxGlobal, fyGlobal;
	PetscReal      **fx, **fy, forceOnProcess[2];

	ierr = DMCompositeGetAccess(lambdaPack, lambda, NULL, &fGlobal); CHKERRQ(ierr);
	ierr = MatMult(ET, fGlobal, regularizedForce);

	ierr = DMCompositeGetAccess(qPack, regularizedForce, &fxGlobal, &fyGlobal); CHKERRQ(ierr);

	// x-direction
	ierr = DMDAVecGetArray(uda, fxGlobal, &fx); CHKERRQ(ierr);
	ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
	forceOnProcess[0] = 0;
	for(PetscInt j=nstart; j<nstart+n; j++)
	{
		for(PetscInt i=mstart; i<mstart+m; i++)
		{
			forceOnProcess[0] += mesh->dy[j] * fx[j][i];
		}
	}
	ierr = DMDAVecRestoreArray(uda, fxGlobal, &fx); CHKERRQ(ierr);

	// y-direction
	ierr = DMDAVecGetArray(vda, fyGlobal, &fy); CHKERRQ(ierr);
	ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
	forceOnProcess[1] = 0;
	for(PetscInt j=nstart; j<nstart+n; j++)
	{
		for(PetscInt i=mstart; i<mstart+m; i++)
		{
			forceOnProcess[1] += mesh->dx[i] * fy[j][i];
		}
	}
	ierr = DMDAVecRestoreArray(vda, fyGlobal, &fy); CHKERRQ(ierr);

	ierr = DMCompositeRestoreAccess(qPack, regularizedForce, &fxGlobal, &fyGlobal); CHKERRQ(ierr);
	ierr = DMCompositeRestoreAccess(lambdaPack, lambda, NULL, &fGlobal); CHKERRQ(ierr);

	ierr = MPI_Reduce(forceOnProcess, force, 2, MPIU_REAL, MPI_SUM, 0, MPI_COMM_WORLD); CHKERRQ(ierr);

	return 0;
}

template <>
PetscErrorCode TairaColoniusSolver<3>::calculateForce()
{
	PetscErrorCode ierr;
	PetscInt       mstart, nstart, pstart, m, n, p;
	Vec            fGlobal, fxGlobal, fyGlobal, fzGlobal;
	PetscReal      ***fx, ***fy, ***fz, forceOnProcess[3];

	ierr = DMCompositeGetAccess(lambdaPack, lambda, NULL, &fGlobal); CHKERRQ(ierr);
	ierr = MatMult(ET, fGlobal, regularizedForce);

	ierr = DMCompositeGetAccess(qPack, regularizedForce, &fxGlobal, &fyGlobal, &fzGlobal); CHKERRQ(ierr);
	
	// x-direction
	ierr = DMDAVecGetArray(uda, fxGlobal, &fx); CHKERRQ(ierr);
	ierr = DMDAGetCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
	forceOnProcess[0] = 0;
	for(PetscInt k=pstart; k<pstart+p; k++)
	{
		for(PetscInt j=nstart; j<nstart+n; j++)
		{
			for(PetscInt i=mstart; i<mstart+m; i++)
			{
				forceOnProcess[0] += mesh->dy[j]*mesh->dz[k] * fx[k][j][i];
			}
		}
	}
	ierr = DMDAVecRestoreArray(uda, fxGlobal, &fx); CHKERRQ(ierr);

	// y-direction
	ierr = DMDAVecGetArray(vda, fyGlobal, &fy); CHKERRQ(ierr);
	ierr = DMDAGetCorners(vda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
	forceOnProcess[1] = 0;
	for(PetscInt k=pstart; k<pstart+p; k++)
	{
		for(PetscInt j=nstart; j<nstart+n; j++)
		{
			for(PetscInt i=mstart; i<mstart+m; i++)
			{
				forceOnProcess[1] += mesh->dz[k]*mesh->dx[i] * fy[k][j][i];
			}
		}
	}
	ierr = DMDAVecRestoreArray(vda, fyGlobal, &fy); CHKERRQ(ierr);

	// z-direction
	ierr = DMDAVecGetArray(wda, fzGlobal, &fz); CHKERRQ(ierr);
	ierr = DMDAGetCorners(wda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
	forceOnProcess[2] = 0;
	for(PetscInt k=pstart; k<pstart+p; k++)
	{
		for(PetscInt j=nstart; j<nstart+n; j++)
		{
			for(PetscInt i=mstart; i<mstart+m; i++)
			{
				forceOnProcess[2] += mesh->dx[i]*mesh->dy[j] * fz[k][j][i];
			}
		}
	}
	ierr = DMDAVecRestoreArray(wda, fzGlobal, &fz); CHKERRQ(ierr);

	ierr = DMCompositeRestoreAccess(qPack, regularizedForce, &fxGlobal, &fyGlobal, &fzGlobal); CHKERRQ(ierr);
	ierr = DMCompositeRestoreAccess(lambdaPack, lambda, NULL, &fGlobal); CHKERRQ(ierr);

	ierr = MPI_Reduce(forceOnProcess, force, 3, MPIU_REAL, MPI_SUM, 0, MPI_COMM_WORLD); CHKERRQ(ierr);

	return 0;
}