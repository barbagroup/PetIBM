template <>
void NavierStokesSolver<2>::createLocalToGlobalMappingsPhi()
{
	PetscErrorCode ierr;
	PetscInt       m, n, i, j, mstart, nstart;
	PetscReal      **lp;
	PetscInt       localIdx;

	// get the range of the vector in the current process
	ierr = VecGetOwnershipRange(phi, &localIdx, NULL); CHKERRV(ierr);

	// populate local vector with the global indices
	// values outside the domain are never accessed and hence not set
	// P
	ierr = DMCreateLocalVector(pda, &pMapping); CHKERRV(ierr);
	ierr = DMDAVecGetArray(pda, pMapping, &lp); CHKERRV(ierr);
	ierr = DMDAGetCorners(pda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRV(ierr);
	for(j=nstart; j<nstart+n; j++)
	{
		for(i=mstart; i<mstart+m; i++)
		{
			lp[j][i] = localIdx;
			localIdx++;
		}
	}
	ierr = DMDAVecRestoreArray(pda, pMapping, &lp); CHKERRV(ierr);

	// scatter from local to local to obtain correct values in ghost cells
	// P
	ierr = DMDALocalToLocalBegin(pda, pMapping, INSERT_VALUES, pMapping); CHKERRV(ierr);
	ierr = DMDALocalToLocalEnd(pda, pMapping, INSERT_VALUES, pMapping); CHKERRV(ierr);
}

template <>
void NavierStokesSolver<3>::createLocalToGlobalMappingsPhi()
{
	PetscErrorCode ierr;
	PetscInt       m, n, p, i, j, k, mstart, nstart, pstart;
	PetscReal      ***lp;
	PetscInt       localIdx;

	// get the range of the vector in the current process
	ierr = VecGetOwnershipRange(phi, &localIdx, NULL); CHKERRV(ierr);

	// populate local vector with the global indices
	// values outside the domain are never accessed and hence not set
	// P
	ierr = DMCreateLocalVector(pda, &pMapping); CHKERRV(ierr);
	ierr = DMDAVecGetArray(pda, pMapping, &lp); CHKERRV(ierr);
	ierr = DMDAGetCorners(pda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRV(ierr);
	for(k=pstart; k<pstart+p; k++)
	{
		for(j=nstart; j<nstart+n; j++)
		{
			for(i=mstart; i<mstart+m; i++)
			{
				lp[k][j][i] = localIdx;
				localIdx++;
			}
		}
	}
	ierr = DMDAVecRestoreArray(pda, pMapping, &lp); CHKERRV(ierr);

	// scatter from local to local to obtain correct values in ghost cells
	// P
	ierr = DMDALocalToLocalBegin(pda, pMapping, INSERT_VALUES, pMapping); CHKERRV(ierr);
	ierr = DMDALocalToLocalEnd(pda, pMapping, INSERT_VALUES, pMapping); CHKERRV(ierr);
}