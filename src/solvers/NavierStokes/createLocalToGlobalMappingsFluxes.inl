template <>
void NavierStokesSolver<2>::createLocalToGlobalMappingsFluxes()
{
	PetscErrorCode ierr;
	PetscInt       m, n, i, j, mstart, nstart;
	PetscReal      **lx, **ly;
	PetscInt       localIdx;

	// get the range of the vector in the current process
	ierr = VecGetOwnershipRange(q, &localIdx, NULL); CHKERRV(ierr);

	// populate local vectors with the global indices
	// set value to -1 if the cell is outside the domain
	// U
	ierr = DMCreateLocalVector(uda, &uMapping); CHKERRV(ierr);
	ierr = DMDAVecGetArray(uda, uMapping, &lx); CHKERRV(ierr);
	ierr = DMDAGetGhostCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRV(ierr);
	for(j=nstart; j<nstart+n; j++)
	{
		for(i=mstart; i<mstart+m; i++)
		{
			lx[j][i] = -1;
			if(i>mstart && i<mstart+m-1 && j>nstart && j<nstart+n-1)
			{
				lx[j][i] = localIdx;
				localIdx++;
			}
		}
	}
	ierr = DMDAVecRestoreArray(uda, uMapping, &lx); CHKERRV(ierr);

	// V
	ierr = DMCreateLocalVector(vda, &vMapping); CHKERRV(ierr);
	ierr = DMDAVecGetArray(vda, vMapping, &ly); CHKERRV(ierr);
	ierr = DMDAGetGhostCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRV(ierr);
	for(j=nstart; j<nstart+n; j++)
	{
		for(i=mstart; i<mstart+m; i++)
		{
			ly[j][i] = -1;
			if(i>mstart && i<mstart+m-1 && j>nstart && j<nstart+n-1)
			{
				ly[j][i] = localIdx;
				localIdx++;
			}	
		}
	}
	ierr = DMDAVecRestoreArray(vda, vMapping, &ly); CHKERRV(ierr);

	// scatter from local to local to obtain correct values in ghost cells
	// U
	ierr = DMDALocalToLocalBegin(uda, uMapping, INSERT_VALUES, uMapping); CHKERRV(ierr);
	ierr = DMDALocalToLocalEnd(uda, uMapping, INSERT_VALUES, uMapping); CHKERRV(ierr);
	// V
	ierr = DMDALocalToLocalBegin(vda, vMapping, INSERT_VALUES, vMapping); CHKERRV(ierr);
	ierr = DMDALocalToLocalEnd(vda, vMapping, INSERT_VALUES, vMapping); CHKERRV(ierr);
}

template <>
void NavierStokesSolver<3>::createLocalToGlobalMappingsFluxes()		
{
	PetscErrorCode ierr;
	PetscInt       i, j, k, m, n, p, mstart, nstart, pstart;
	PetscReal      ***lx, ***ly, ***lz;
	PetscInt       localIdx;

	// get the range of the vector in the current process
	ierr = VecGetOwnershipRange(q, &localIdx, NULL); CHKERRV(ierr);

	// populate local vectors with the global indices
	// set value to -1 if the cell is outside the domain
	// U
	ierr = DMCreateLocalVector(uda, &uMapping); CHKERRV(ierr);
	ierr = DMDAVecGetArray(uda, uMapping, &lx); CHKERRV(ierr);
	ierr = DMDAGetGhostCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRV(ierr);
	for(k=pstart; k<pstart+p; k++)
	{
		for(j=nstart; j<nstart+n; j++)
		{
			for(i=mstart; i<mstart+m; i++)
			{
				lx[k][j][i] = -1;
				if(i>mstart && i<mstart+m-1 && j>nstart && j<nstart+n-1 && k>pstart && k<pstart+p-1)
				{
					lx[k][j][i] = localIdx;
					localIdx++;
				}
			}
		}
	}
	ierr = DMDAVecRestoreArray(uda, uMapping, &lx); CHKERRV(ierr);
	// V
	ierr = DMCreateLocalVector(vda, &vMapping); CHKERRV(ierr);
	ierr = DMDAVecGetArray(vda, vMapping, &ly); CHKERRV(ierr);
	ierr = DMDAGetGhostCorners(vda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRV(ierr);
	for(k=pstart; k<pstart+p; k++)
	{
		for(j=nstart; j<nstart+n; j++)
		{
			for(i=mstart; i<mstart+m; i++)
			{
				ly[k][j][i] = -1;
				if(i>mstart && i<mstart+m-1 && j>nstart && j<nstart+n-1 && k>pstart && k<pstart+p-1)
				{
					ly[k][j][i] = localIdx;
					localIdx++;
				}
			}
		}
	}
	ierr = DMDAVecRestoreArray(vda, vMapping, &ly); CHKERRV(ierr);
	// W
	ierr = DMCreateLocalVector(wda, &wMapping); CHKERRV(ierr);
	ierr = DMDAVecGetArray(wda, wMapping, &lz); CHKERRV(ierr);
	ierr = DMDAGetGhostCorners(wda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRV(ierr);
	for(k=pstart; k<pstart+p; k++)
	{
		for(j=nstart; j<nstart+n; j++)
		{
			for(i=mstart; i<mstart+m; i++)
			{
				lz[k][j][i] = -1;
				if(i>mstart && i<mstart+m-1 && j>nstart && j<nstart+n-1 && k>pstart && k<pstart+p-1)
				{
					lz[k][j][i] = localIdx;
					localIdx++;
				}
			}
		}
	}
	ierr = DMDAVecRestoreArray(wda, wMapping, &lz); CHKERRV(ierr);

	// scatter from local to local to obtain correct values in ghost cells
	// U
	ierr = DMDALocalToLocalBegin(uda, uMapping, INSERT_VALUES, uMapping); CHKERRV(ierr);
	ierr = DMDALocalToLocalEnd(uda, uMapping, INSERT_VALUES, uMapping); CHKERRV(ierr);
	// V
	ierr = DMDALocalToLocalBegin(vda, vMapping, INSERT_VALUES, vMapping); CHKERRV(ierr);
	ierr = DMDALocalToLocalEnd(vda, vMapping, INSERT_VALUES, vMapping); CHKERRV(ierr);
	// W
	ierr = DMDALocalToLocalBegin(wda, wMapping, INSERT_VALUES, wMapping); CHKERRV(ierr);
	ierr = DMDALocalToLocalEnd(wda, wMapping, INSERT_VALUES, wMapping); CHKERRV(ierr);
}