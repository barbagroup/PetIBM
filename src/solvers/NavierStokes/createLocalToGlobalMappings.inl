template <>
void NavierStokesSolver<2>::createLocalToGlobalMappings()
{
	PetscErrorCode ierr;
	PetscInt       m, n, i, j, mstart, nstart;
	Vec            globalX, globalY;
	PetscReal      **gx, **gy, **lx, **ly;
	PetscInt       localIdx;

	ierr = DMCreateGlobalVector(pack, &globalIndices); CHKERRV(ierr);
	ierr = DMCreateLocalVector(uda, &uMapping); CHKERRV(ierr);
	ierr = DMCreateLocalVector(vda, &vMapping); CHKERRV(ierr);

	// initialise global vector
	ierr = VecGetOwnershipRange(globalIndices, &localIdx, NULL); CHKERRV(ierr);
	ierr = DMCompositeGetAccess(pack, globalIndices, &globalX, &globalY); CHKERRV(ierr);

	// U
	ierr = DMDAVecGetArray(uda, globalX, &gx); CHKERRV(ierr);
	ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRV(ierr);
	for(j=nstart; j<nstart+n; j++)
	{
		for(i=mstart; i<mstart+m; i++)
		{
			gx[j][i] = localIdx;
			localIdx++;
		}
	}
	ierr = DMDAVecRestoreArray(uda, globalX, &gx); CHKERRV(ierr);

	// V
	ierr = DMDAVecGetArray(vda, globalY, &gy); CHKERRV(ierr);
	ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRV(ierr);
	for(j=nstart; j<nstart+n; j++)
	{
		for(i=mstart; i<mstart+m; i++)
		{
			gy[j][i] = localIdx;
			localIdx++;
		}
	}
	ierr = DMDAVecRestoreArray(vda, globalY, &gy); CHKERRV(ierr);

	// initialise local vector
	// set values to -1
	// U
	ierr = DMDAVecGetArray(uda, uMapping, &lx); CHKERRV(ierr);
	ierr = DMDAGetGhostCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRV(ierr);
	for(j=nstart; j<nstart+n; j++)
	{
		for(i=mstart; i<mstart+m; i++)
		{
			lx[j][i] = -1;
		}
	}
	ierr = DMDAVecRestoreArray(uda, uMapping, &lx); CHKERRV(ierr);

	// V
	ierr = DMDAVecGetArray(vda, vMapping, &ly); CHKERRV(ierr);
	ierr = DMDAGetGhostCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRV(ierr);
	for(j=nstart; j<nstart+n; j++)
	{
		for(i=mstart; i<mstart+m; i++)
		{
			ly[j][i] = -1;
		}
	}
	ierr = DMDAVecRestoreArray(vda, vMapping, &ly); CHKERRV(ierr);

	// scatter from global to local vectors
	ierr = DMGlobalToLocalBegin(uda, globalX, INSERT_VALUES, uMapping); CHKERRV(ierr);
	ierr = DMGlobalToLocalEnd(uda, globalX, INSERT_VALUES, uMapping); CHKERRV(ierr);

	ierr = DMGlobalToLocalBegin(vda, globalY, INSERT_VALUES, vMapping); CHKERRV(ierr);
	ierr = DMGlobalToLocalEnd(vda, globalY, INSERT_VALUES, vMapping); CHKERRV(ierr);

	ierr = DMCompositeRestoreAccess(pack, globalIndices, &globalX, &globalY); CHKERRV(ierr);

	ierr = VecDestroy(&globalIndices);
}

template <>
void NavierStokesSolver<3>::createLocalToGlobalMappings()		
{
	PetscErrorCode ierr;
	PetscInt       i, j, k, m, n, p, mstart, nstart, pstart;
	Vec            globalX, globalY, globalZ;
	PetscReal      ***gx, ***gy, ***gz, ***lx, ***ly, ***lz;
	PetscInt       localIdx;

	ierr = DMCreateGlobalVector(pack, &globalIndices); CHKERRV(ierr);
	ierr = DMCreateLocalVector(uda, &uMapping); CHKERRV(ierr);
	ierr = DMCreateLocalVector(vda, &vMapping); CHKERRV(ierr);
	ierr = DMCreateLocalVector(wda, &wMapping); CHKERRV(ierr);

	// initialise global vector
	ierr = VecGetOwnershipRange(globalIndices, &localIdx, NULL); CHKERRV(ierr);
	ierr = DMCompositeGetAccess(pack, globalIndices, &globalX, &globalY, &globalZ); CHKERRV(ierr);
	// U
	ierr = DMDAVecGetArray(uda, globalX, &gx); CHKERRV(ierr);
	ierr = DMDAGetCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRV(ierr);
	for(k=pstart; k<pstart+p; k++)
	{
		for(j=nstart; j<nstart+n; j++)
		{
			for(i=mstart; i<mstart+m; i++)
			{
				gx[k][j][i] = localIdx;
				localIdx++;
			}
		}
	}
	ierr = DMDAVecRestoreArray(uda, globalX, &gx); CHKERRV(ierr);
	// V
	ierr = DMDAVecGetArray(vda, globalY, &gy); CHKERRV(ierr);
	ierr = DMDAGetCorners(vda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRV(ierr);
	for(k=pstart; k<pstart+p; k++)
	{
		for(j=nstart; j<nstart+n; j++)
		{
			for(i=mstart; i<mstart+m; i++)
			{
				gy[k][j][i] = localIdx;
				localIdx++;
			}
		}
	}
	ierr = DMDAVecRestoreArray(vda, globalY, &gy); CHKERRV(ierr);
	// W
	ierr = DMDAVecGetArray(wda, globalZ, &gz); CHKERRV(ierr);
	ierr = DMDAGetCorners(wda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRV(ierr);
	for(k=pstart; k<pstart+p; k++)
	{
		for(j=nstart; j<nstart+n; j++)
		{
			for(i=mstart; i<mstart+m; i++)
			{
				gz[k][j][i] = localIdx;
				localIdx++;
			}
		}
	}
	ierr = DMDAVecRestoreArray(wda, globalZ, &gz); CHKERRV(ierr);

	// initialise local vector
	// set values to -1
	// U
	ierr = DMDAVecGetArray(uda, uMapping, &lx); CHKERRV(ierr);
	ierr = DMDAGetGhostCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRV(ierr);
	for(k=pstart; k<pstart+p; k++)
	{
		for(j=nstart; j<nstart+n; j++)
		{
			for(i=mstart; i<mstart+m; i++)
			{
				lx[k][j][i] = -1;
			}
		}
	}
	ierr = DMDAVecRestoreArray(uda, uMapping, &lx); CHKERRV(ierr);
	// V
	ierr = DMDAVecGetArray(vda, vMapping, &ly); CHKERRV(ierr);
	ierr = DMDAGetGhostCorners(vda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRV(ierr);
	for(k=pstart; k<pstart+p; k++)
	{
		for(j=nstart; j<nstart+n; j++)
		{
			for(i=mstart; i<mstart+m; i++)
			{
				ly[k][j][i] = -1;
			}
		}
	}
	ierr = DMDAVecRestoreArray(vda, vMapping, &ly); CHKERRV(ierr);
	// W
	ierr = DMDAVecGetArray(wda, wMapping, &lz); CHKERRV(ierr);
	ierr = DMDAGetGhostCorners(wda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRV(ierr);
	for(k=pstart; k<pstart+p; k++)
	{
		for(j=nstart; j<nstart+n; j++)
		{
			for(i=mstart; i<mstart+m; i++)
			{
				lz[k][j][i] = -1;
			}
		}
	}
	ierr = DMDAVecRestoreArray(wda, wMapping, &lz); CHKERRV(ierr);

	// scatter from global to local vectors
	// U
	ierr = DMGlobalToLocalBegin(uda, globalX, INSERT_VALUES, uMapping); CHKERRV(ierr);
	ierr = DMGlobalToLocalEnd(uda, globalX, INSERT_VALUES, uMapping); CHKERRV(ierr);
	// V
	ierr = DMGlobalToLocalBegin(vda, globalY, INSERT_VALUES, vMapping); CHKERRV(ierr);
	ierr = DMGlobalToLocalEnd(vda, globalY, INSERT_VALUES, vMapping); CHKERRV(ierr);
	// W
	ierr = DMGlobalToLocalBegin(wda, globalZ, INSERT_VALUES, wMapping); CHKERRV(ierr);
	ierr = DMGlobalToLocalEnd(wda, globalZ, INSERT_VALUES, wMapping); CHKERRV(ierr);

	ierr = DMCompositeRestoreAccess(pack, globalIndices, &globalX, &globalY); CHKERRV(ierr);
	ierr = VecDestroy(&globalIndices);
}