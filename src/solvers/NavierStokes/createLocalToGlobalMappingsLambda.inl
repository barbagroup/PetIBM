template <>
PetscErrorCode NavierStokesSolver<2>::createLocalToGlobalMappingsLambda()
{
	PetscErrorCode ierr;
	PetscInt       m, n, i, j, mstart, nstart;
	PetscReal      **lp;
	PetscInt       localIdx;

	// get the range of the vector in the current process
	ierr = VecGetOwnershipRange(lambda, &localIdx, NULL); CHKERRQ(ierr);

	// populate local vector with the global indices
	// values outside the domain are never accessed and hence not set
	// P
	ierr = DMCreateLocalVector(pda, &pMapping); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(pda, pMapping, &lp); CHKERRQ(ierr);
	ierr = DMDAGetCorners(pda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
	for(j=nstart; j<nstart+n; j++)
	{
		for(i=mstart; i<mstart+m; i++)
		{
			lp[j][i] = localIdx;
			localIdx++;
		}
	}
	ierr = DMDAVecRestoreArray(pda, pMapping, &lp); CHKERRQ(ierr);

	// scatter from local to local to obtain correct values in ghost cells
	// P
	ierr = DMDALocalToLocalBegin(pda, pMapping, INSERT_VALUES, pMapping); CHKERRQ(ierr);
	ierr = DMDALocalToLocalEnd(pda, pMapping, INSERT_VALUES, pMapping); CHKERRQ(ierr);

	return 0;
}

template <>
PetscErrorCode NavierStokesSolver<3>::createLocalToGlobalMappingsLambda()
{
	PetscErrorCode ierr;
	PetscInt       m, n, p, i, j, k, mstart, nstart, pstart;
	PetscReal      ***lp;
	PetscInt       localIdx;

	// get the range of the vector in the current process
	ierr = VecGetOwnershipRange(lambda, &localIdx, NULL); CHKERRQ(ierr);

	// populate local vector with the global indices
	// values outside the domain are never accessed and hence not set
	// P
	ierr = DMCreateLocalVector(pda, &pMapping); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(pda, pMapping, &lp); CHKERRQ(ierr);
	ierr = DMDAGetCorners(pda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
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
	ierr = DMDAVecRestoreArray(pda, pMapping, &lp); CHKERRQ(ierr);

	// scatter from local to local to obtain correct values in ghost cells
	// P
	ierr = DMDALocalToLocalBegin(pda, pMapping, INSERT_VALUES, pMapping); CHKERRQ(ierr);
	ierr = DMDALocalToLocalEnd(pda, pMapping, INSERT_VALUES, pMapping); CHKERRQ(ierr);

	return 0;
}