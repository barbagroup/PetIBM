/*! Implementation of the method `createGlobalMappingBodies` of the class `TairaColoniusSolver`.
 * \file createGlobalMappingBodies.inl
 */


/*!
 * \brief Maps local to global indices.
 */
template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::createGlobalMappingBodies()
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  PetscMPIInt numProcs;
  ierr = MPI_Comm_size(PETSC_COMM_WORLD, &numProcs); CHKERRQ(ierr);

  globalIndexMapping.resize(numLagrangianPoints);

  PetscInt globalIndex = 0;
  PetscInt offset;
  for (PetscInt procIdx=0; procIdx<numProcs; procIdx++)
  {
    globalIndex += localNumPhiPoints[procIdx];
    offset = 0;
    for (PetscInt bIdx=0; bIdx<numBodies; bIdx++)
    {
      for (auto i=bodies[bIdx].localIndexPoints[procIdx].begin(); i!=bodies[bIdx].localIndexPoints[procIdx].end(); i++)
      {
        globalIndexMapping[offset + *i] = globalIndex;
        globalIndex += dim;
      }
      offset += bodies[bIdx].localIndexPoints[procIdx].size();
    }
  }

  PetscFunctionReturn(0);
} // createGlobalMappingBodies
