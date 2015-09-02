/***************************************************************************//**
 * \file generateBodyInfo.inl
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the method `generateBodyInfo` of the class `TairaColoniusSolver`.
 */


/**
 * \brief
 */
template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::generateBodyInfo()
{
  return 0;
} // generateBodyInfo


// two-dimensional specialization
template <>
PetscErrorCode TairaColoniusSolver<2>::generateBodyInfo()
{
  PetscErrorCode ierr;
  
  PetscInt numProcs;
  ierr = MPI_Comm_size(PETSC_COMM_WORLD, &numProcs); CHKERRQ(ierr);

  boundaryPointIndices.resize(numProcs);
  numBoundaryPointsOnProcess.resize(numProcs);
  numPhiOnProcess.resize(numProcs);

  globalIndexMapping.resize(bodies[0].numPoints); // why here?

  const PetscInt *plx, *ply;
  ierr = DMDAGetOwnershipRanges(pda, &plx, &ply, NULL); CHKERRQ(ierr);
  PetscInt m, n; // number of procs in each dimension
  ierr = DMDAGetInfo(pda, NULL, NULL, NULL, NULL, &m, &n, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);

  PetscInt xStart, yStart, xEnd, yEnd,
           procIdx;

  yStart = 0;
  for (PetscInt j=0; j<n; j++)
  {
    yEnd = yStart + ply[j];
    xStart = 0;
    for (PetscInt i=0; i<m; i++)
    {
      procIdx = j*m + i;
      xEnd = xStart + plx[i];
      numPhiOnProcess[procIdx] = plx[i]*ply[j];
      for (PetscInt l=0; l<bodies[0].numPoints; l++)
      {
        if (bodies[0].X[l] >= mesh->x[xStart] && bodies[0].X[l] < mesh->x[xEnd] 
            && bodies[0].Y[l] >= mesh->y[yStart] && bodies[0].Y[l] < mesh->y[yEnd])
        {
          numBoundaryPointsOnProcess[procIdx]++;
        }
      }
      boundaryPointIndices[procIdx].reserve(numBoundaryPointsOnProcess[procIdx]);
      for (PetscInt l=0; l<bodies[0].numPoints; l++)
      {
        if (bodies[0].X[l] >= mesh->x[xStart] && bodies[0].X[l] < mesh->x[xEnd] 
            && bodies[0].Y[l] >= mesh->y[yStart] && bodies[0].Y[l] < mesh->y[yEnd])
        {
          boundaryPointIndices[procIdx].push_back(l);
        }
      }
      xStart = xEnd;
    }
    yStart = yEnd;
  }

  return 0;
} // generateBodyInfo


// three-dimensional specialization
template <>
PetscErrorCode TairaColoniusSolver<3>::generateBodyInfo()
{
  PetscErrorCode ierr;

  PetscInt numProcs;  
  ierr = MPI_Comm_size(PETSC_COMM_WORLD, &numProcs); CHKERRQ(ierr);

  boundaryPointIndices.resize(numProcs);
  numBoundaryPointsOnProcess.resize(numProcs);
  numPhiOnProcess.resize(numProcs);

  globalIndexMapping.resize(bodies[0].numPoints);

  const PetscInt *plx, *ply, *plz;
  ierr = DMDAGetOwnershipRanges(pda, &plx, &ply, &plz); CHKERRQ(ierr);
  PetscInt m, n, p; // local number of nodes along each direction
  ierr = DMDAGetInfo(pda, NULL, NULL, NULL, NULL, &m, &n, &p, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);

  PetscInt xStart, yStart, zStart, xEnd, yEnd, zEnd,
           procIdx = 0;
  
  zStart = 0;
  for (PetscInt k=0; k<p; k++)
  {
    zEnd = zStart + plz[k];
    yStart = 0;
    for (PetscInt j=0; j<n; j++)
    {
      yEnd = yStart + ply[j];
      xStart = 0;
      for (PetscInt i=0; i<m; i++)
      {
        procIdx = k*m*n + j*m + i;
        xEnd = xStart + plx[i];
        numPhiOnProcess[procIdx] = plx[i]*ply[j]*plz[k];
        for (PetscInt l=0; l<bodies[0].numPoints; l++)
        {
          if (bodies[0].X[l] >= mesh->x[xStart] && bodies[0].X[l] < mesh->x[xEnd] 
              && bodies[0].Y[l] >= mesh->y[yStart] && bodies[0].Y[l] < mesh->y[yEnd] 
              && bodies[0].Z[l] >= mesh->z[zStart] && bodies[0].Z[l] < mesh->z[zEnd])
          {
            numBoundaryPointsOnProcess[procIdx]++;
          }
        }
        boundaryPointIndices[procIdx].reserve(numBoundaryPointsOnProcess[procIdx]);
        for (PetscInt l=0; l<bodies[0].numPoints; l++)
        {
          if (bodies[0].X[l] >= mesh->x[xStart] && bodies[0].X[l] < mesh->x[xEnd] 
              && bodies[0].Y[l] >= mesh->y[yStart] && bodies[0].Y[l] < mesh->y[yEnd] 
              && bodies[0].Z[l] >= mesh->z[zStart] && bodies[0].Z[l] < mesh->z[zEnd])
          {
            boundaryPointIndices[procIdx].push_back(l);
          }
        }
        xStart = xEnd;
      }
      yStart = yEnd;
    }
    zStart = zEnd;
  }

  return 0;
} // generateBodyInfo