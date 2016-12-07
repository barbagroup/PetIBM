/***************************************************************************//**
 * \file createGlobalMappingBodies.inl
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the method `createGlobalMappingBodies` 
 *        of the class `TairaColoniusSolver`.
 */


/**
 * \brief Maps local to global indices.
 */
template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::createGlobalMappingBodies()
{
  return 0;
} // createGlobalMappingBodies


// two-dimensional specialization
template <>
PetscErrorCode TairaColoniusSolver<2>::createGlobalMappingBodies()
{
  PetscErrorCode ierr;

  PetscMPIInt numProcs;
  ierr = MPI_Comm_size(PETSC_COMM_WORLD, &numProcs); CHKERRQ(ierr);

  PetscInt globalIndex = 0;
  for (PetscInt procIdx=0; procIdx<numProcs; procIdx++)
  {
    globalIndex += numPhiOnProcess[procIdx];
    for (auto i=boundaryPointIndices[procIdx].begin(); i!=boundaryPointIndices[procIdx].end(); i++)
    {
      globalIndexMapping[*i] = globalIndex;
      globalIndex+=2;
    }
  }

  return 0;
} // createGlobalMappingBodies


// three-dimensional specialization
template <>
PetscErrorCode TairaColoniusSolver<3>::createGlobalMappingBodies()
{
  PetscErrorCode ierr;

  PetscMPIInt numProcs;
  ierr = MPI_Comm_size(PETSC_COMM_WORLD, &numProcs); CHKERRQ(ierr);

  PetscInt globalIndex = 0;
  for (PetscInt procIdx=0; procIdx<numProcs; procIdx++)
  {
    globalIndex += numPhiOnProcess[procIdx];
    for (auto i=boundaryPointIndices[procIdx].begin(); i!=boundaryPointIndices[procIdx].end(); i++)
    {
      globalIndexMapping[*i] = globalIndex;
      globalIndex+=3;
    }
  }

  return 0;
} // createGlobalMappingBodies