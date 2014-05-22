template <>
void NavierStokesSolver<2>::generateBNQ()
{
	PetscErrorCode ierr;
	PetscInt       i, j;
	PetscInt       mstart, nstart, m, n;
	PetscInt       *d_nnz, *o_nnz;
	PetscInt       qStart, qEnd, phiStart, phiEnd, qLocalSize, phiLocalSize;
	PetscInt       localIdx;
	PetscReal      **pGlobalIdx;
	PetscInt       row, cols[2];
	PetscReal      values[2] = {-1.0, 1.0};
	
	ierr = VecGetOwnershipRange(q, &qStart, &qEnd); CHKERRV(ierr);
	qLocalSize = qEnd-qStart;

	// create arrays to store nnz values
	ierr = PetscMalloc(qLocalSize*sizeof(PetscInt), &d_nnz); CHKERRV(ierr);
	ierr = PetscMalloc(qLocalSize*sizeof(PetscInt), &o_nnz); CHKERRV(ierr);
	
	// ownership range of PHI
	ierr = VecGetOwnershipRange(phi, &phiStart, &phiEnd); CHKERRV(ierr);
	phiLocalSize = phiEnd-phiStart;

	ierr = DMDAVecGetArray(pda, pMapping, &pGlobalIdx); CHKERRV(ierr);

	// determine the number of non-zeros in each row
	// in the diagonal and off-diagonal portions of the matrix
	localIdx = 0;
	// U
	ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRV(ierr);
	for(j=nstart; j<nstart+n; j++)
	{
		for(i=mstart; i<mstart+m; i++)
		{
			cols[0] = pGlobalIdx[j][i];
			cols[1] = pGlobalIdx[j][i+1];
			countNumNonZeros(cols, 2, phiStart, phiEnd, d_nnz[localIdx], o_nnz[localIdx]);
			localIdx++;
		}
	}
	// V
	ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRV(ierr);
	for(j=nstart; j<nstart+n; j++)
	{
		for(i=mstart; i<mstart+m; i++)
		{
			cols[0] = pGlobalIdx[j][i];
			cols[1] = pGlobalIdx[j+1][i];
			countNumNonZeros(cols, 2, phiStart, phiEnd, d_nnz[localIdx], o_nnz[localIdx]);
			localIdx++;
		}
	}
	
	// allocate memory for the matrix
	ierr = MatCreate(PETSC_COMM_WORLD, &BNQ); CHKERRV(ierr);
	ierr = MatSetType(BNQ, MATMPIAIJ); CHKERRV(ierr);
	ierr = MatSetSizes(BNQ, qLocalSize, phiLocalSize, PETSC_DETERMINE, PETSC_DETERMINE); CHKERRV(ierr);
	ierr = MatMPIAIJSetPreallocation(BNQ, 0, d_nnz, 0, o_nnz); CHKERRV(ierr);

	// deallocate d_nnz and o_nnz
	ierr = PetscFree(d_nnz); CHKERRV(ierr);
	ierr = PetscFree(o_nnz); CHKERRV(ierr);

	// assemble matrix Q
	localIdx = 0;
	// U
	ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRV(ierr);
	for(j=nstart; j<nstart+n; j++)
	{
		for(i=mstart; i<mstart+m; i++)
		{
			row = localIdx + qStart;
			cols[0] = pGlobalIdx[j][i];
			cols[1] = pGlobalIdx[j][i+1];
			ierr = MatSetValues(BNQ, 1, &row, 2, cols, values, INSERT_VALUES); CHKERRV(ierr);
			localIdx++;
		}
	}
	// V
	ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRV(ierr);
	for(j=nstart; j<nstart+n; j++)
	{
		for(i=mstart; i<mstart+m; i++)
		{
			row = localIdx + qStart;
			cols[0] = pGlobalIdx[j][i];
			cols[1] = pGlobalIdx[j+1][i];
			ierr = MatSetValues(BNQ, 1, &row, 2, cols, values, INSERT_VALUES); CHKERRV(ierr);
			localIdx++;
		}
	}
	ierr = DMDAVecRestoreArray(pda, pMapping, &pGlobalIdx); CHKERRV(ierr);

	ierr = MatAssemblyBegin(BNQ, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
	ierr = MatAssemblyEnd(BNQ, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);

	ierr = MatTranspose(BNQ, MAT_INITIAL_MATRIX, &QT); CHKERRV(ierr);
	ierr = MatDiagonalScale(BNQ, BN, NULL); CHKERRV(ierr);
}

template <>
void NavierStokesSolver<3>::generateBNQ()
{
	PetscErrorCode ierr;
	PetscInt       i, j, k;
	PetscInt       mstart, nstart, pstart, m, n, p;
	PetscInt       *d_nnz, *o_nnz;
	PetscInt       qStart, qEnd, phiStart, phiEnd, qLocalSize, phiLocalSize;
	PetscInt       localIdx;
	PetscReal      ***pGlobalIdx;
	PetscInt       row, cols[2];
	PetscReal      values[2] = {-1.0, 1.0};	
	
	ierr = VecGetOwnershipRange(q, &qStart, &qEnd); CHKERRV(ierr);
	qLocalSize = qEnd-qStart;

	// create arrays to store nnz values
	ierr = PetscMalloc(qLocalSize*sizeof(PetscInt), &d_nnz); CHKERRV(ierr);
	ierr = PetscMalloc(qLocalSize*sizeof(PetscInt), &o_nnz); CHKERRV(ierr);
	
	// ownership range of PHI
	ierr = VecGetOwnershipRange(phi, &phiStart, &phiEnd); CHKERRV(ierr);
	phiLocalSize = phiEnd-phiStart;

	ierr = DMDAVecGetArray(pda, pMapping, &pGlobalIdx); CHKERRV(ierr);

	// determine the number of non-zeros in each row
	// in the diagonal and off-diagonal portions of the matrix
	localIdx = 0;
	// U
	ierr = DMDAGetCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRV(ierr);
	for(k=pstart; k<pstart+p; k++)
	{
		for(j=nstart; j<nstart+n; j++)
		{
			for(i=mstart; i<mstart+m; i++)
			{
				cols[0] = pGlobalIdx[k][j][i];
				cols[1] = pGlobalIdx[k][j][i+1];
				countNumNonZeros(cols, 2, phiStart, phiEnd, d_nnz[localIdx], o_nnz[localIdx]);
				localIdx++;
			}
		}
	}
	// V
	ierr = DMDAGetCorners(vda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRV(ierr);
	for(k=pstart; k<pstart+p; k++)
	{
		for(j=nstart; j<nstart+n; j++)
		{
			for(i=mstart; i<mstart+m; i++)
			{
				cols[0] = pGlobalIdx[k][j][i];
				cols[1] = pGlobalIdx[k][j+1][i];
				countNumNonZeros(cols, 2, phiStart, phiEnd, d_nnz[localIdx], o_nnz[localIdx]);
				localIdx++;
			}
		}
	}
	// W
	ierr = DMDAGetCorners(wda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRV(ierr);
	for(k=pstart; k<pstart+p; k++)
	{
		for(j=nstart; j<nstart+n; j++)
		{
			for(i=mstart; i<mstart+m; i++)
			{
				cols[0] = pGlobalIdx[k][j][i];
				cols[1] = pGlobalIdx[k+1][j][i];
				countNumNonZeros(cols, 2, phiStart, phiEnd, d_nnz[localIdx], o_nnz[localIdx]);
				localIdx++;
			}
		}
	}
	
	// allocate memory for the matrix
	ierr = MatCreate(PETSC_COMM_WORLD, &BNQ); CHKERRV(ierr);
	ierr = MatSetType(BNQ, MATMPIAIJ); CHKERRV(ierr);
	ierr = MatSetSizes(BNQ, qLocalSize, phiLocalSize, PETSC_DETERMINE, PETSC_DETERMINE); CHKERRV(ierr);
	ierr = MatMPIAIJSetPreallocation(BNQ, 0, d_nnz, 0, o_nnz); CHKERRV(ierr);

	// deallocate d_nnz and o_nnz
	ierr = PetscFree(d_nnz); CHKERRV(ierr);
	ierr = PetscFree(o_nnz); CHKERRV(ierr);

	// assemble matrix Q
	localIdx = 0;
	// U
	ierr = DMDAGetCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRV(ierr);
	for(k=pstart; k<pstart+p; k++)
	{
		for(j=nstart; j<nstart+n; j++)
		{
			for(i=mstart; i<mstart+m; i++)
			{
				row = localIdx + qStart;
				cols[0] = pGlobalIdx[k][j][i];
				cols[1] = pGlobalIdx[k][j][i+1];
				ierr = MatSetValues(BNQ, 1, &row, 2, cols, values, INSERT_VALUES); CHKERRV(ierr);
				localIdx++;
			}
		}
	}
	// V
	ierr = DMDAGetCorners(vda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRV(ierr);
	for(k=pstart; k<pstart+p; k++)
	{
		for(j=nstart; j<nstart+n; j++)
		{
			for(i=mstart; i<mstart+m; i++)
			{
				row = localIdx + qStart;
				cols[0] = pGlobalIdx[k][j][i];
				cols[1] = pGlobalIdx[k][j+1][i];
				ierr = MatSetValues(BNQ, 1, &row, 2, cols, values, INSERT_VALUES); CHKERRV(ierr);
				localIdx++;
			}
		}
	}
	// W
	ierr = DMDAGetCorners(wda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRV(ierr);
	for(k=pstart; k<pstart+p; k++)
	{
		for(j=nstart; j<nstart+n; j++)
		{
			for(i=mstart; i<mstart+m; i++)
			{
				row = localIdx + qStart;
				cols[0] = pGlobalIdx[k][j][i];
				cols[1] = pGlobalIdx[k+1][j][i];
				ierr = MatSetValues(BNQ, 1, &row, 2, cols, values, INSERT_VALUES); CHKERRV(ierr);
				localIdx++;
			}
		}
	}
	ierr = DMDAVecRestoreArray(pda, pMapping, &pGlobalIdx); CHKERRV(ierr);

	ierr = MatAssemblyBegin(BNQ, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
	ierr = MatAssemblyEnd(BNQ, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);

	ierr = MatTranspose(BNQ, MAT_INITIAL_MATRIX, &QT); CHKERRV(ierr);
	ierr = MatDiagonalScale(BNQ, BN, NULL); CHKERRV(ierr);
}