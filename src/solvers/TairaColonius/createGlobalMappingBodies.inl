template <>
void TairaColoniusSolver<2>::createGlobalMappingBodies()
{
	PetscErrorCode ierr;
	PetscInt       m, n;
	PetscInt       lambdaStart, lambdaEnd;
	PetscInt       numProcs;

	MPI_Comm_size(PETSC_COMM_WORLD, &numProcs);

	ierr = VecGetOwnershipRange(lambda, &lambdaStart, &lambdaEnd); CHKERRV(ierr);
	ierr = DMDAGetCorners(pda, NULL, NULL, NULL, &m, &n, NULL); CHKERRV(ierr);
	startGlobalIndex = lambdaStart + m*n;

	MPI_Barrier(PETSC_COMM_WORLD);
	MPI_Allgather(&startGlobalIndex, 1, MPIU_INT, &startGlobalIndices.front(), 1, MPIU_INT, PETSC_COMM_WORLD);

	PetscInt globalIndex;

	for(PetscInt j=0; j<numProcs; j++)
	{
		globalIndex = startGlobalIndices[j];
		for(auto i=boundaryPointIndices[j].begin(); i!=boundaryPointIndices[j].end(); i++)
		{
			bodyGlobalIndices[*i] = globalIndex;
			globalIndex++;
		}
	}
}

template <>
void TairaColoniusSolver<3>::createGlobalMappingBodies()
{
	PetscErrorCode ierr;
	PetscInt       m, n, p;
	PetscInt       lambdaStart, lambdaEnd;
	PetscInt       numProcs;

	MPI_Comm_size(PETSC_COMM_WORLD, &numProcs);

	ierr = VecGetOwnershipRange(lambda, &lambdaStart, &lambdaEnd); CHKERRV(ierr);
	ierr = DMDAGetCorners(pda, NULL, NULL, NULL, &m, &n, &p); CHKERRV(ierr);
	startGlobalIndex = lambdaStart + m*n*p;

	MPI_Barrier(PETSC_COMM_WORLD);
	MPI_Allgather(&startGlobalIndex, 1, MPIU_INT, &startGlobalIndices.front(), 1, MPIU_INT, PETSC_COMM_WORLD);

	PetscInt globalIndex;

	for(PetscInt j=0; j<numProcs; j++)
	{
		globalIndex = startGlobalIndices[j];
		for(auto i=boundaryPointIndices[j].begin(); i!=boundaryPointIndices[j].end(); i++)
		{
			bodyGlobalIndices[*i] = globalIndex;
			globalIndex++;
		}
	}
}