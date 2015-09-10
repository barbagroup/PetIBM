/***************************************************************************//**
 * \file createLocalToGlobalMappingsLambda.inl
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the method `createLocalToGlobalMappingsLambda` 
 *        of the class `NavierStokesSolver`.
 */


/**
 * \brief Maps local multi-dimensional indices to global index for lambda.
 *
 * Vectors stored as distributed arrays can be accessed using multi-dimensional 
 * arrays on every process, with each index referring to the numbering along 
 * each cartesian direction. The elements of the vector also have a global
 * ordering. This function generates the map from the local multi-dimensional
 * indexing pf the pressure variables to the global indices of the vector 
 * \f$ \lambda \f$.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::createLocalToGlobalMappingsLambda()
{
  return 0;
} // createLocalToGlobalMappingsLambda


// two-dimensional specialization
template <>
PetscErrorCode NavierStokesSolver<2>::createLocalToGlobalMappingsLambda()
{
  PetscErrorCode ierr;
  PetscInt i, j,           // loop indices
           m, n,           // local number of nodes along each direction
           mstart, nstart; // starting indices

  // get global index of first local element of pressure
  PetscInt globalIdx;
  ierr = VecGetOwnershipRange(lambda, &globalIdx, NULL); CHKERRQ(ierr);

  // populate local vector with global indices
  // values outside the domain are never accessed and hence not set
  ierr = DMCreateLocalVector(pda, &pMapping); CHKERRQ(ierr);
  PetscReal **pMappingArray;
  ierr = DMDAVecGetArray(pda, pMapping, &pMappingArray); CHKERRQ(ierr);
  ierr = DMDAGetCorners(pda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
  for (j=nstart; j<nstart+n; j++)
  {
    for (i=mstart; i<mstart+m; i++)
    {
      pMappingArray[j][i] = globalIdx;
      globalIdx++;
    }
  }
  ierr = DMDAVecRestoreArray(pda, pMapping, &pMappingArray); CHKERRQ(ierr);

  // scatter from local to local to obtain correct values in ghost cells
  ierr = DMLocalToLocalBegin(pda, pMapping, INSERT_VALUES, pMapping); CHKERRQ(ierr);
  ierr = DMLocalToLocalEnd(pda, pMapping, INSERT_VALUES, pMapping); CHKERRQ(ierr);

  return 0;
} // createLocalToGlobalMappingsLambda


// three-dimensional specialization
template <>
PetscErrorCode NavierStokesSolver<3>::createLocalToGlobalMappingsLambda()
{
  PetscErrorCode ierr;
  PetscInt i, j, k,                // loop indices
           m, n, p,                // local number of nodes along each direction
           mstart, nstart, pstart; // starting indices

  // get global index of first local element of pressure
  PetscInt globalIdx;
  ierr = VecGetOwnershipRange(lambda, &globalIdx, NULL); CHKERRQ(ierr);

  // populate local vector with the global indices
  // values outside the domain are never accessed and hence not set
  ierr = DMCreateLocalVector(pda, &pMapping); CHKERRQ(ierr);
  PetscReal ***pMappingArray;
  ierr = DMDAVecGetArray(pda, pMapping, &pMappingArray); CHKERRQ(ierr);
  ierr = DMDAGetCorners(pda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  for (k=pstart; k<pstart+p; k++)
  {
    for (j=nstart; j<nstart+n; j++)
    {
      for (i=mstart; i<mstart+m; i++)
      {
        pMappingArray[k][j][i] = globalIdx;
        globalIdx++;
      }
    }
  }
  ierr = DMDAVecRestoreArray(pda, pMapping, &pMappingArray); CHKERRQ(ierr);

  // scatter from local to local to obtain correct values in ghost cells
  ierr = DMLocalToLocalBegin(pda, pMapping, INSERT_VALUES, pMapping); CHKERRQ(ierr);
  ierr = DMLocalToLocalEnd(pda, pMapping, INSERT_VALUES, pMapping); CHKERRQ(ierr);

  return 0;
} // createLocalToGlobalMappingsLambda