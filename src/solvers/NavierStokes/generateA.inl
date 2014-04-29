template <>
void NavierStokesSolver<2>::generateA()
{
	PetscErrorCode ierr;
	PetscInt       i, j, l;
	PetscInt       mstart, nstart;
	PetscInt       m, n;
	PetscInt       *d_nnz, *o_nnz;
	PetscInt       qStart, qEnd, qLocalSize, pStart, pEnd, pLocalSize;
	PetscInt       row, col, cols[5], localCols[5];
	PetscInt       localIdx;
	PetscReal      x, y, value, values[2];
	ISLocalToGlobalMapping *ltogs;

	// map local sub-DM (including ghost) indices to packed global indices
	ierr = DMCompositeGetISLocalToGlobalMappings(pack, &ltogs); CHKERRV(ierr);

	// ownership range of flux vector
	ierr = VecGetOwnershipRange(q, &qStart, &qEnd); CHKERRV(ierr);
	qLocalSize = qEnd-qStart;

	ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRV(ierr);
	ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Number of local U: %d\n", m*n); CHKERRV(ierr);
	for(j=0; j<n; j++)
	{
		for(i=0; i<m; i++)
		{
			localIdx = (j+1)*(m+2) + (i+1); // this assumes the presence of ghost cells

			localCols[0] = localIdx-(m+2);
			localCols[1] = localIdx-1;
			localCols[2] = localIdx;
			localCols[3] = localIdx+1;
			localCols[4] = localIdx+(m+2);

			ierr = ISLocalToGlobalMappingApply(ltogs[0], 5, localCols, cols); CHKERRV(ierr);
			ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "%3d |%3d:%3d |%3d:%3d |%3d:%3d |%3d:%3d |%3d:%3d\n", cols[2], localCols[0], cols[0], localCols[1], cols[1], localCols[2], cols[2], localCols[3], cols[3], localCols[4], cols[4]);
		}
	}

	ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRV(ierr);
	ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Number of local V: %d\n", m*n); CHKERRV(ierr);
	for(j=0; j<n; j++)
	{
		for(i=0; i<m; i++)
		{
			localIdx = (j+1)*(m+2) + (i+1);

			localCols[0] = localIdx-(m+2);
			localCols[1] = localIdx-1;
			localCols[2] = localIdx;
			localCols[3] = localIdx+1;
			localCols[4] = localIdx+(m+2);

			ierr = ISLocalToGlobalMappingApply(ltogs[1], 5, localCols, cols); CHKERRV(ierr);
			ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "%3d |%3d:%3d | %3d:%3d |%3d:%3d |%3d:%3d |%3d:%3d\n", cols[2], localCols[0], cols[0], localCols[1], cols[1], localCols[2], cols[2], localCols[3], cols[3], localCols[4], cols[4]);
		}
	}
	ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Total U and V    : %d\n\n", qLocalSize); CHKERRV(ierr);
	ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD);

	ierr = ISLocalToGlobalMappingDestroy(&ltogs[0]); CHKERRV(ierr);
	ierr = ISLocalToGlobalMappingDestroy(&ltogs[1]); CHKERRV(ierr);
	ierr = PetscFree(ltogs); CHKERRV(ierr);
}

template <>
void NavierStokesSolver<3>::generateA()
{
	
}