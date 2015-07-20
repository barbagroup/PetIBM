template <>
PetscErrorCode TairaColoniusSolver<2>::createGlobalMappingBodies()
{
	PetscErrorCode ierr;
	PetscInt       numProcs;
	PetscInt       globalIndex;

	ierr = MPI_Comm_size(PETSC_COMM_WORLD, &numProcs); CHKERRQ(ierr);

	globalIndex = 0;
	for(PetscInt procIdx=0; procIdx<numProcs; procIdx++)
	{
		globalIndex += numPhiOnProcess[procIdx];
		for(auto i=boundaryPointIndices[procIdx].begin(); i!=boundaryPointIndices[procIdx].end(); i++)
		{
			globalIndexMapping[*i] = globalIndex;
			globalIndex+=2;
		}
	}

	return 0;
} // createGlobalMappingBodies

template <>
PetscErrorCode TairaColoniusSolver<3>::createGlobalMappingBodies()
{
	PetscErrorCode ierr;
	PetscInt       numProcs;
	PetscInt       globalIndex;

	ierr = MPI_Comm_size(PETSC_COMM_WORLD, &numProcs); CHKERRQ(ierr);

	globalIndex = 0;
	for(PetscInt procIdx=0; procIdx<numProcs; procIdx++)
	{
		globalIndex += numPhiOnProcess[procIdx];
		for(auto i=boundaryPointIndices[procIdx].begin(); i!=boundaryPointIndices[procIdx].end(); i++)
		{
			globalIndexMapping[*i] = globalIndex;
			globalIndex+=3;
		}
	}

	return 0;
} // createGlobalMappingBodies