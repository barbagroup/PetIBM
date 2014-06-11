template <>
PetscErrorCode TairaColoniusSolver<2>::createGlobalMappingBodies()
{
	PetscInt       numProcs;
	PetscInt       globalIndex;

	MPI_Comm_size(PETSC_COMM_WORLD, &numProcs);

	globalIndex = 0;
	for(PetscInt procIdx=0; procIdx<numProcs; procIdx++)
	{
		globalIndex += numPhiOnProcess[procIdx];
		for(auto i=boundaryPointIndices[procIdx].begin(); i!=boundaryPointIndices[procIdx].end(); i++)
		{
			globalIndexMapping[*i] = globalIndex;
			globalIndex++;
		}
		globalIndex += boundaryPointIndices[procIdx].size();
	}

	return 0;
}

template <>
PetscErrorCode TairaColoniusSolver<3>::createGlobalMappingBodies()
{
	PetscInt       numProcs;
	PetscInt       globalIndex;

	MPI_Comm_size(PETSC_COMM_WORLD, &numProcs);

	globalIndex = 0;
	for(PetscInt procIdx=0; procIdx<numProcs; procIdx++)
	{
		globalIndex += numPhiOnProcess[procIdx];
		for(auto i=boundaryPointIndices[procIdx].begin(); i!=boundaryPointIndices[procIdx].end(); i++)
		{
			globalIndexMapping[*i] = globalIndex;
			globalIndex++;
		}
		globalIndex += 2*boundaryPointIndices[procIdx].size();
	}

	return 0;
}