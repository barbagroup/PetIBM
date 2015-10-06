/***************************************************************************//**
 * \file generateBNQ.inl
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the method `generateBNQ` of the class `NavierStokesSolver`.
 */


/**
 * \brief Assembles the matrix \f$ B^N Q \f$.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::generateBNQ()
{
  return 0;
} // generateBNQ


// two-dimensional specialization
template <>
PetscErrorCode NavierStokesSolver<2>::generateBNQ()
{
  PetscErrorCode ierr;

  PetscInt i, j,           // loop indices 
           m, n,           // local number of nodes along each direction
           mstart, nstart; // startting indices

  PetscInt localIdx;
  PetscInt row, cols[2];
  PetscReal values[2] = {-1.0, 1.0};
  
  // get local size of fluxes vector
  PetscInt qStart, qEnd, qLocalSize;
  ierr = VecGetOwnershipRange(q, &qStart, &qEnd); CHKERRQ(ierr);
  qLocalSize = qEnd-qStart;

  // create arrays to store number of non-zeros (nnz) values
  PetscInt *d_nnz, // nnz on diagonal
           *o_nnz; // nnz off diagonal
  ierr = PetscMalloc(qLocalSize*sizeof(PetscInt), &d_nnz); CHKERRQ(ierr);
  ierr = PetscMalloc(qLocalSize*sizeof(PetscInt), &o_nnz); CHKERRQ(ierr);
  
  // get local size of lambda vector
  PetscInt lambdaStart, lambdaEnd, lambdaLocalSize;
  ierr = VecGetOwnershipRange(lambda, &lambdaStart, &lambdaEnd); CHKERRQ(ierr);
  lambdaLocalSize = lambdaEnd-lambdaStart;

  PetscReal **pMappingArray;
  ierr = DMDAVecGetArray(pda, pMapping, &pMappingArray); CHKERRQ(ierr);

  // determine nnz in matrix BNQ row by row
  localIdx = 0;
  // rows corresponding to fluxes in x-direction
  ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
  for (j=nstart; j<nstart+n; j++)
  {
    for (i=mstart; i<mstart+m; i++)
    {
      cols[0] = pMappingArray[j][i];
      cols[1] = pMappingArray[j][i+1];
      countNumNonZeros(cols, 2, lambdaStart, lambdaEnd, d_nnz[localIdx], o_nnz[localIdx]);
      localIdx++;
    }
  }
  // rows corresponding to fluxes in y-direction
  ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
  for (j=nstart; j<nstart+n; j++)
  {
    for (i=mstart; i<mstart+m; i++)
    {
      cols[0] = pMappingArray[j][i];
      cols[1] = pMappingArray[j+1][i];
      countNumNonZeros(cols, 2, lambdaStart, lambdaEnd, d_nnz[localIdx], o_nnz[localIdx]);
      localIdx++;
    }
  }
  
  // allocate memory for matrix BNQ
  ierr = MatCreate(PETSC_COMM_WORLD, &BNQ); CHKERRQ(ierr);
  ierr = MatSetSizes(BNQ, qLocalSize, lambdaLocalSize, PETSC_DETERMINE, PETSC_DETERMINE); CHKERRQ(ierr);
  ierr = MatSetFromOptions(BNQ); CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(BNQ, 0, d_nnz); CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(BNQ, 0, d_nnz, 0, o_nnz); CHKERRQ(ierr);
  ierr = MatSetOption(BNQ, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE); CHKERRQ(ierr);

  // deallocate d_nnz and o_nnz
  ierr = PetscFree(d_nnz); CHKERRQ(ierr);
  ierr = PetscFree(o_nnz); CHKERRQ(ierr);

  // assemble matrix Q row by row
  localIdx = 0;
  // rows corresponding to fluxes in x-direction
  ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
  for (j=nstart; j<nstart+n; j++)
  {
    for (i=mstart; i<mstart+m; i++)
    {
      row = localIdx + qStart;
      cols[0] = pMappingArray[j][i];
      cols[1] = pMappingArray[j][i+1];
      ierr = MatSetValues(BNQ, 1, &row, 2, cols, values, INSERT_VALUES); CHKERRQ(ierr);
      localIdx++;
    }
  }
  // rows corresponding to fluxes in y-direction
  ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
  for (j=nstart; j<nstart+n; j++)
  {
    for (i=mstart; i<mstart+m; i++)
    {
      row = localIdx + qStart;
      cols[0] = pMappingArray[j][i];
      cols[1] = pMappingArray[j+1][i];
      ierr = MatSetValues(BNQ, 1, &row, 2, cols, values, INSERT_VALUES); CHKERRQ(ierr);
      localIdx++;
    }
  }
  ierr = DMDAVecRestoreArray(pda, pMapping, &pMappingArray); CHKERRQ(ierr);

  ierr = MatAssemblyBegin(BNQ, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(BNQ, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  // compute matrix QT
  ierr = MatTranspose(BNQ, MAT_INITIAL_MATRIX, &QT); CHKERRQ(ierr);
  // scale Q to get BNQ
  ierr = MatDiagonalScale(BNQ, BN, NULL); CHKERRQ(ierr);

  return 0;
} // generateBNQ


// three-dimensional specialization
template <>
PetscErrorCode NavierStokesSolver<3>::generateBNQ()
{
  PetscErrorCode ierr;

  PetscInt i, j, k,                // loop indices
           m, n, p,                // local number of nodes along each direction
           mstart, nstart, pstart; // starting indices
  
  PetscInt localIdx;
  PetscInt row, cols[2];
  PetscReal values[2] = {-1.0, 1.0}; 
  
  // get local size of fluxes vector
  PetscInt qStart, qEnd, qLocalSize;
  ierr = VecGetOwnershipRange(q, &qStart, &qEnd); CHKERRQ(ierr);
  qLocalSize = qEnd-qStart;

  // create arrays to store number of non-zeros (nnz) values
  PetscInt *d_nnz, // nnz on diagonal
           *o_nnz; // nnz off diagonal
  ierr = PetscMalloc(qLocalSize*sizeof(PetscInt), &d_nnz); CHKERRQ(ierr);
  ierr = PetscMalloc(qLocalSize*sizeof(PetscInt), &o_nnz); CHKERRQ(ierr);
  
  // get local size of lambda vector
  PetscInt lambdaStart, lambdaEnd, lambdaLocalSize;
  ierr = VecGetOwnershipRange(lambda, &lambdaStart, &lambdaEnd); CHKERRQ(ierr);
  lambdaLocalSize = lambdaEnd-lambdaStart;

  PetscReal ***pMappingArray;
  ierr = DMDAVecGetArray(pda, pMapping, &pMappingArray); CHKERRQ(ierr);

  // determine nnz in BNQ row by row
  localIdx = 0;
  // rows corresponding to fluxes in x-direction
  ierr = DMDAGetCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  for (k=pstart; k<pstart+p; k++)
  {
    for (j=nstart; j<nstart+n; j++)
    {
      for (i=mstart; i<mstart+m; i++)
      {
        cols[0] = pMappingArray[k][j][i];
        cols[1] = pMappingArray[k][j][i+1];
        countNumNonZeros(cols, 2, lambdaStart, lambdaEnd, d_nnz[localIdx], o_nnz[localIdx]);
        localIdx++;
      }
    }
  }
  // rows corresponding to fluxes in y-direction
  ierr = DMDAGetCorners(vda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  for (k=pstart; k<pstart+p; k++)
  {
    for (j=nstart; j<nstart+n; j++)
    {
      for (i=mstart; i<mstart+m; i++)
      {
        cols[0] = pMappingArray[k][j][i];
        cols[1] = pMappingArray[k][j+1][i];
        countNumNonZeros(cols, 2, lambdaStart, lambdaEnd, d_nnz[localIdx], o_nnz[localIdx]);
        localIdx++;
      }
    }
  }
  // rows corresponding to fluxes in z-direction
  ierr = DMDAGetCorners(wda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  for (k=pstart; k<pstart+p; k++)
  {
    for (j=nstart; j<nstart+n; j++)
    {
      for (i=mstart; i<mstart+m; i++)
      {
        cols[0] = pMappingArray[k][j][i];
        cols[1] = pMappingArray[k+1][j][i];
        countNumNonZeros(cols, 2, lambdaStart, lambdaEnd, d_nnz[localIdx], o_nnz[localIdx]);
        localIdx++;
      }
    }
  }
  
  // allocate memory for matrix BNQ
  ierr = MatCreate(PETSC_COMM_WORLD, &BNQ); CHKERRQ(ierr);
  ierr = MatSetSizes(BNQ, qLocalSize, lambdaLocalSize, PETSC_DETERMINE, PETSC_DETERMINE); CHKERRQ(ierr);
  ierr = MatSetFromOptions(BNQ); CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(BNQ, 0, d_nnz); CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(BNQ, 0, d_nnz, 0, o_nnz); CHKERRQ(ierr);
  ierr = MatSetOption(BNQ, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE); CHKERRQ(ierr);

  // deallocate d_nnz and o_nnz
  ierr = PetscFree(d_nnz); CHKERRQ(ierr);
  ierr = PetscFree(o_nnz); CHKERRQ(ierr);

  // assemble matrix Q row by row
  localIdx = 0;
  // rows corresponding to fluxes in x-direction
  ierr = DMDAGetCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  for (k=pstart; k<pstart+p; k++)
  {
    for (j=nstart; j<nstart+n; j++)
    {
      for (i=mstart; i<mstart+m; i++)
      {
        row = localIdx + qStart;
        cols[0] = pMappingArray[k][j][i];
        cols[1] = pMappingArray[k][j][i+1];
        ierr = MatSetValues(BNQ, 1, &row, 2, cols, values, INSERT_VALUES); CHKERRQ(ierr);
        localIdx++;
      }
    }
  }
  // rows corresponding to fluxes in y-direction
  ierr = DMDAGetCorners(vda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  for (k=pstart; k<pstart+p; k++)
  {
    for (j=nstart; j<nstart+n; j++)
    {
      for (i=mstart; i<mstart+m; i++)
      {
        row = localIdx + qStart;
        cols[0] = pMappingArray[k][j][i];
        cols[1] = pMappingArray[k][j+1][i];
        ierr = MatSetValues(BNQ, 1, &row, 2, cols, values, INSERT_VALUES); CHKERRQ(ierr);
        localIdx++;
      }
    }
  }
  // rows corresponding to fluxes in z-direction
  ierr = DMDAGetCorners(wda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  for (k=pstart; k<pstart+p; k++)
  {
    for (j=nstart; j<nstart+n; j++)
    {
      for (i=mstart; i<mstart+m; i++)
      {
        row = localIdx + qStart;
        cols[0] = pMappingArray[k][j][i];
        cols[1] = pMappingArray[k+1][j][i];
        ierr = MatSetValues(BNQ, 1, &row, 2, cols, values, INSERT_VALUES); CHKERRQ(ierr);
        localIdx++;
      }
    }
  }
  ierr = DMDAVecRestoreArray(pda, pMapping, &pMappingArray); CHKERRQ(ierr);

  ierr = MatAssemblyBegin(BNQ, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(BNQ, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  // compute matrix QT
  ierr = MatTranspose(BNQ, MAT_INITIAL_MATRIX, &QT); CHKERRQ(ierr);
  // scale Q to get BNQ
  ierr = MatDiagonalScale(BNQ, BN, NULL); CHKERRQ(ierr);

  return 0;
} // generateBNQ
