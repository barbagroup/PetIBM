/***************************************************************************//**
* Vectors stored as distributed arrays can be accessed using multi-dimensional 
* arrays on every process, with each index referring to the numbering along 
* each cartesian direction. The elements of the vector also have a global
* ordering. This function generates the map from the multi-dimensional
* indexing of each of the local flux vectors `qx`, `qy` and `qz`, to the global 
* indices of the composite flux vector `q`.
*/
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::createLocalToGlobalMappingsFluxes()
{
	return 0;
}

template <>
PetscErrorCode NavierStokesSolver<2>::createLocalToGlobalMappingsFluxes()
{
	PetscErrorCode ierr;
	PetscInt       m, n, i, j, mstart, nstart;
	PetscReal      **lx, **ly;
	PetscInt       localIdx;

	// get the range of the vector in the current process
	ierr = VecGetOwnershipRange(q, &localIdx, NULL); CHKERRQ(ierr);

	// populate local vectors with the global indices
	// set value to -1 if the cell is outside the domain
	// U
	ierr = DMCreateLocalVector(uda, &uMapping); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(uda, uMapping, &lx); CHKERRQ(ierr);
	ierr = DMDAGetGhostCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
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
	ierr = DMDAVecRestoreArray(uda, uMapping, &lx); CHKERRQ(ierr);

	// V
	ierr = DMCreateLocalVector(vda, &vMapping); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(vda, vMapping, &ly); CHKERRQ(ierr);
	ierr = DMDAGetGhostCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
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
	ierr = DMDAVecRestoreArray(vda, vMapping, &ly); CHKERRQ(ierr);

	// scatter from local to local to obtain correct values in ghost cells
	// U
	ierr = DMDALocalToLocalBegin(uda, uMapping, INSERT_VALUES, uMapping); CHKERRQ(ierr);
	ierr = DMDALocalToLocalEnd(uda, uMapping, INSERT_VALUES, uMapping); CHKERRQ(ierr);
	// V
	ierr = DMDALocalToLocalBegin(vda, vMapping, INSERT_VALUES, vMapping); CHKERRQ(ierr);
	ierr = DMDALocalToLocalEnd(vda, vMapping, INSERT_VALUES, vMapping); CHKERRQ(ierr);

	return 0;
}

template <>
PetscErrorCode NavierStokesSolver<3>::createLocalToGlobalMappingsFluxes()		
{
	PetscErrorCode ierr;
	PetscInt       i, j, k, m, n, p, mstart, nstart, pstart;
	PetscReal      ***lx, ***ly, ***lz;
	PetscInt       localIdx;

	// get the range of the vector in the current process
	ierr = VecGetOwnershipRange(q, &localIdx, NULL); CHKERRQ(ierr);

	// populate local vectors with the global indices
	// set value to -1 if the cell is outside the domain
	// U
	ierr = DMCreateLocalVector(uda, &uMapping); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(uda, uMapping, &lx); CHKERRQ(ierr);
	ierr = DMDAGetGhostCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
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
	ierr = DMDAVecRestoreArray(uda, uMapping, &lx); CHKERRQ(ierr);
	// V
	ierr = DMCreateLocalVector(vda, &vMapping); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(vda, vMapping, &ly); CHKERRQ(ierr);
	ierr = DMDAGetGhostCorners(vda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
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
	ierr = DMDAVecRestoreArray(vda, vMapping, &ly); CHKERRQ(ierr);
	// W
	ierr = DMCreateLocalVector(wda, &wMapping); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(wda, wMapping, &lz); CHKERRQ(ierr);
	ierr = DMDAGetGhostCorners(wda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
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
	ierr = DMDAVecRestoreArray(wda, wMapping, &lz); CHKERRQ(ierr);

	// scatter from local to local to obtain correct values in ghost cells
	// U
	ierr = DMDALocalToLocalBegin(uda, uMapping, INSERT_VALUES, uMapping); CHKERRQ(ierr);
	ierr = DMDALocalToLocalEnd(uda, uMapping, INSERT_VALUES, uMapping); CHKERRQ(ierr);
	// V
	ierr = DMDALocalToLocalBegin(vda, vMapping, INSERT_VALUES, vMapping); CHKERRQ(ierr);
	ierr = DMDALocalToLocalEnd(vda, vMapping, INSERT_VALUES, vMapping); CHKERRQ(ierr);
	// W
	ierr = DMDALocalToLocalBegin(wda, wMapping, INSERT_VALUES, wMapping); CHKERRQ(ierr);
	ierr = DMDALocalToLocalEnd(wda, wMapping, INSERT_VALUES, wMapping); CHKERRQ(ierr);

	return 0;
}