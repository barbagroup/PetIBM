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
  for (PetscInt procIdx=0; procIdx<numProcs; procIdx++)
  {
    globalIndex += localNumPhiPoints[procIdx];
    for (PetscInt bIdx=0; bIdx<numBodies; bIdx++)
    {
      for (auto i=bodies[bIdx].localIndexPoints[procIdx].begin(); i!=bodies[bIdx].localIndexPoints[procIdx].end(); i++)
      {
        globalIndexMapping[*i] = globalIndex;
        globalIndex += dim;
      }
    }
  }

  PetscFunctionReturn(0);
} // createGlobalMappingBodies
