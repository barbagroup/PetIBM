/*! Implementation of the method `generateBodyInfo` of the class `TairaColoniusSolver`.
 * \file getBodyInfo.inl
 */


/*!
 * \brief Counts the number of Lagrangian points on each process.
 */
template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::getBodyInfo(Body<dim> body,
                                                     std::vector<PetscInt> &localNumberPoints,
                                                     std::vector<std::vector<PetscInt> > &localIndexPoints)
{
  return 0;
} // getBodyInfo


// two-dimensional specialization
template <>
PetscErrorCode TairaColoniusSolver<2>::getBodyInfo(Body<2> body,
                                                   std::vector<PetscInt> &localNumberPoints,
                                                   std::vector<std::vector<PetscInt> > &localIndexPoints)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  
  PetscMPIInt numProcs;
  ierr = MPI_Comm_size(PETSC_COMM_WORLD, &numProcs); CHKERRQ(ierr);

  localNumberPoints.resize(numProcs);
  localIndexPoints.resize(numProcs);
  numPhiOnProcess.resize(numProcs);
  globalIndexMapping.resize(body.numPoints); // why here?  

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
      PetscReal box[4] = {mesh->x[xStart], mesh->x[xEnd],
                          mesh->y[yStart], mesh->y[yEnd]};
      ierr = body.getIndexPointsInBox(box, localIndexPoints[procIdx]); CHKERRQ(ierr);
      localNumberPoints[procIdx] = localIndexPoints[procIdx].size();
      xStart = xEnd;
    }
    yStart = yEnd;
  }

  PetscFunctionReturn(0);
} // getBodyInfo


// three-dimensional specialization
template <>
PetscErrorCode TairaColoniusSolver<3>::getBodyInfo(Body<3> body,
                                                   std::vector<PetscInt> &localNumberPoints,
                                                   std::vector<std::vector<PetscInt> > &localIndexPoints)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  PetscMPIInt numProcs;  
  ierr = MPI_Comm_size(PETSC_COMM_WORLD, &numProcs); CHKERRQ(ierr);

  localNumberPoints.resize(numProcs);
  localIndexPoints.resize(numProcs);
  numPhiOnProcess.resize(numProcs);
  globalIndexMapping.resize(body.numPoints); // why here?  

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
        PetscReal box[6] = {mesh->x[xStart], mesh->x[xEnd],
                            mesh->y[yStart], mesh->y[yEnd],
                            mesh->z[zStart], mesh->z[zEnd]};
        // ierr = bodies[0].getIndexPointsInBox(box, &boundaryPointIndices[procIdx]); CHKERRQ(ierr);
        // numBoundaryPointsOnProcess[procIdx] = boundaryPointIndices[procIdx].size();
        ierr = body.getIndexPointsInBox(box, localIndexPoints[procIdx]); CHKERRQ(ierr);
        localNumberPoints[procIdx] = localIndexPoints[procIdx].size();
        xStart = xEnd;
      }
      yStart = yEnd;
    }
    zStart = zEnd;
  }

  PetscFunctionReturn(0);
} // getBodyInfo
