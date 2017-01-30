/*! Implementation of the method `setLocalNumPhiPoints` of the class `TairaColoniusSolver`.
 * \file setLocalNumPhiPoints.inl
 */


/*!
 * \brief Calculates and sets the number of pressure points on each process.
 */
template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::setLocalNumPhiPoints()
{
  return 0;
} // setLocalNumPhiPoints


// two-dimensional specialization
template <>
PetscErrorCode TairaColoniusSolver<2>::setLocalNumPhiPoints()
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  PetscMPIInt numProcs;
  ierr = MPI_Comm_size(PETSC_COMM_WORLD, &numProcs); CHKERRQ(ierr);

  localNumPhiPoints.resize(numProcs);

  const PetscInt *plx, *ply;
  ierr = DMDAGetOwnershipRanges(pda, &plx, &ply, NULL); CHKERRQ(ierr);
  PetscInt m, n; // number of procs in each dimension
  ierr = DMDAGetInfo(pda, NULL, NULL, NULL, NULL, &m, &n, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);

  PetscInt procIdx;
  for (PetscInt j=0; j<n; j++)
  {
    for (PetscInt i=0; i<m; i++)
    {
      procIdx = j*m + i;
      localNumPhiPoints[procIdx] = plx[i]*ply[j];
    }
  }

  PetscFunctionReturn(0);
} // setLocalNumPhiPoints


// three-dimensional specialization
template <>
PetscErrorCode TairaColoniusSolver<3>::setLocalNumPhiPoints()
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  PetscMPIInt numProcs;  
  ierr = MPI_Comm_size(PETSC_COMM_WORLD, &numProcs); CHKERRQ(ierr);

  localNumPhiPoints.resize(numProcs);

  const PetscInt *plx, *ply, *plz;
  ierr = DMDAGetOwnershipRanges(pda, &plx, &ply, &plz); CHKERRQ(ierr);
  PetscInt m, n, p; // local number of nodes along each direction
  ierr = DMDAGetInfo(pda, NULL, NULL, NULL, NULL, &m, &n, &p, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);

  PetscInt procIdx;
  for (PetscInt k=0; k<p; k++)
  {
    for (PetscInt j=0; j<n; j++)
    {
      for (PetscInt i=0; i<m; i++)
      {
        procIdx = k*m*n + j*m + i;
        localNumPhiPoints[procIdx] = plx[i]*ply[j]*plz[k];
      }
    }
  }

  PetscFunctionReturn(0);
} // setLocalNumPhiPoints
