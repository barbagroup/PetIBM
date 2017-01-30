/*! Implementation of the method `setLocalIndexPointsBodies` of the class `TairaColoniusSolver`.
 * \file setLocalIndexPointsBodies.inl
 */


/*!
 * \brief Calculates and sets the number of Lagrangian points on each process.
 */
template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::setLocalIndexPointsBodies()
{
  return 0;
} // setLocalIndexPointsBodies


// two-dimensional specialization
template <>
PetscErrorCode TairaColoniusSolver<2>::setLocalIndexPointsBodies()
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

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
      PetscReal box[4] = {mesh->x[xStart], mesh->x[xEnd],
                          mesh->y[yStart], mesh->y[yEnd]};
      for (PetscInt l=0; l<numBodies; l++)
      {
        ierr = bodies[l].setLocalIndexPoints(procIdx, box); CHKERRQ(ierr);
      }
      xStart = xEnd;
    }
    yStart = yEnd;
  }

  PetscFunctionReturn(0);
} // setLocalIndexPointsBodies


// three-dimensional specialization
template <>
PetscErrorCode TairaColoniusSolver<3>::setLocalIndexPointsBodies()
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

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
        PetscReal box[6] = {mesh->x[xStart], mesh->x[xEnd],
                            mesh->y[yStart], mesh->y[yEnd],
                            mesh->z[zStart], mesh->z[zEnd]};
        for (PetscInt l=0; l<numBodies; l++)
        {
          ierr = bodies[l].setLocalIndexPoints(procIdx, box); CHKERRQ(ierr);
        }
        xStart = xEnd;
      }
      yStart = yEnd;
    }
    zStart = zEnd;
  }

  PetscFunctionReturn(0);
} // setLocalIndexPointsBodies
