/***************************************************************************//**
 * \file generateBNQ.inl
 * \author Anush Krishnan (anush@bu.edu)
 * \brief
 */


/**
 * \brief
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
  PetscInt i, j, 
           m, n, 
           mstart, nstart;

  PetscInt localIdx;
  PetscInt row, cols[2];
  PetscReal values[2] = {-1.0, 1.0};
  
  PetscInt qStart, qEnd, qLocalSize;
  ierr = VecGetOwnershipRange(q, &qStart, &qEnd); CHKERRQ(ierr);
  qLocalSize = qEnd-qStart;

  // create arrays to store nnz values
  PetscInt *d_nnz, *o_nnz;
  ierr = PetscMalloc(qLocalSize*sizeof(PetscInt), &d_nnz); CHKERRQ(ierr);
  ierr = PetscMalloc(qLocalSize*sizeof(PetscInt), &o_nnz); CHKERRQ(ierr);
  
  // ownership range of lambda
  PetscInt lambdaStart, lambdaEnd, lambdaLocalSize;
  ierr = VecGetOwnershipRange(lambda, &lambdaStart, &lambdaEnd); CHKERRQ(ierr);
  lambdaLocalSize = lambdaEnd-lambdaStart;

  PetscReal **pGlobalIdx;
  ierr = DMDAVecGetArray(pda, pMapping, &pGlobalIdx); CHKERRQ(ierr);

  // determine the number of non-zeros in each row
  // in the diagonal and off-diagonal portions of the matrix
  localIdx = 0;
  // non-zeros for x-momentum equation
  ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
  for (j=nstart; j<nstart+n; j++)
  {
    for (i=mstart; i<mstart+m; i++)
    {
      cols[0] = pGlobalIdx[j][i];
      cols[1] = pGlobalIdx[j][i+1];
      countNumNonZeros(cols, 2, lambdaStart, lambdaEnd, d_nnz[localIdx], o_nnz[localIdx]);
      localIdx++;
    }
  }
  // non-zeros for y-momentum equation
  ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
  for (j=nstart; j<nstart+n; j++)
  {
    for (i=mstart; i<mstart+m; i++)
    {
      cols[0] = pGlobalIdx[j][i];
      cols[1] = pGlobalIdx[j+1][i];
      countNumNonZeros(cols, 2, lambdaStart, lambdaEnd, d_nnz[localIdx], o_nnz[localIdx]);
      localIdx++;
    }
  }
  
  // allocate memory for the matrix
  ierr = MatCreate(PETSC_COMM_WORLD, &BNQ); CHKERRQ(ierr);
  ierr = MatSetSizes(BNQ, qLocalSize, lambdaLocalSize, PETSC_DETERMINE, PETSC_DETERMINE); CHKERRQ(ierr);
  ierr = MatSetFromOptions(BNQ); CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(BNQ, 0, d_nnz); CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(BNQ, 0, d_nnz, 0, o_nnz); CHKERRQ(ierr);

  // deallocate d_nnz and o_nnz
  ierr = PetscFree(d_nnz); CHKERRQ(ierr);
  ierr = PetscFree(o_nnz); CHKERRQ(ierr);

  // assemble matrix Q
  localIdx = 0;
  // x-momentum coefficients
  ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
  for (j=nstart; j<nstart+n; j++)
  {
    for (i=mstart; i<mstart+m; i++)
    {
      row = localIdx + qStart;
      cols[0] = pGlobalIdx[j][i];
      cols[1] = pGlobalIdx[j][i+1];
      ierr = MatSetValues(BNQ, 1, &row, 2, cols, values, INSERT_VALUES); CHKERRQ(ierr);
      localIdx++;
    }
  }
  // y-momentum coefficients
  ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
  for (j=nstart; j<nstart+n; j++)
  {
    for (i=mstart; i<mstart+m; i++)
    {
      row = localIdx + qStart;
      cols[0] = pGlobalIdx[j][i];
      cols[1] = pGlobalIdx[j+1][i];
      ierr = MatSetValues(BNQ, 1, &row, 2, cols, values, INSERT_VALUES); CHKERRQ(ierr);
      localIdx++;
    }
  }
  ierr = DMDAVecRestoreArray(pda, pMapping, &pGlobalIdx); CHKERRQ(ierr);

  ierr = MatAssemblyBegin(BNQ, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(BNQ, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  ierr = MatTranspose(BNQ, MAT_INITIAL_MATRIX, &QT); CHKERRQ(ierr);
  ierr = MatDiagonalScale(BNQ, BN, NULL); CHKERRQ(ierr);

  return 0;
} // generateBNQ


// three-dimensional specialization
template <>
PetscErrorCode NavierStokesSolver<3>::generateBNQ()
{
  PetscErrorCode ierr;
  PetscInt i, j, k, 
           m, n, p, 
           mstart, nstart, pstart;
  
  PetscInt localIdx;
  PetscInt row, cols[2];
  PetscReal values[2] = {-1.0, 1.0}; 
  
  PetscInt qStart, qEnd, qLocalSize;
  ierr = VecGetOwnershipRange(q, &qStart, &qEnd); CHKERRQ(ierr);
  qLocalSize = qEnd-qStart;

  // create arrays to store nnz values
  PetscInt *d_nnz, *o_nnz;
  ierr = PetscMalloc(qLocalSize*sizeof(PetscInt), &d_nnz); CHKERRQ(ierr);
  ierr = PetscMalloc(qLocalSize*sizeof(PetscInt), &o_nnz); CHKERRQ(ierr);
  
  // ownership range of lambda
  PetscInt lambdaStart, lambdaEnd, lambdaLocalSize;
  ierr = VecGetOwnershipRange(lambda, &lambdaStart, &lambdaEnd); CHKERRQ(ierr);
  lambdaLocalSize = lambdaEnd-lambdaStart;

  PetscReal ***pGlobalIdx;
  ierr = DMDAVecGetArray(pda, pMapping, &pGlobalIdx); CHKERRQ(ierr);

  // determine the number of non-zeros in each row
  // in the diagonal and off-diagonal portions of the matrix
  localIdx = 0;
  // non-zeros for x-momentum equation
  ierr = DMDAGetCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  for (k=pstart; k<pstart+p; k++)
  {
    for (j=nstart; j<nstart+n; j++)
    {
      for (i=mstart; i<mstart+m; i++)
      {
        cols[0] = pGlobalIdx[k][j][i];
        cols[1] = pGlobalIdx[k][j][i+1];
        countNumNonZeros(cols, 2, lambdaStart, lambdaEnd, d_nnz[localIdx], o_nnz[localIdx]);
        localIdx++;
      }
    }
  }
  // non-zeros for y-momentum equation
  ierr = DMDAGetCorners(vda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  for (k=pstart; k<pstart+p; k++)
  {
    for (j=nstart; j<nstart+n; j++)
    {
      for (i=mstart; i<mstart+m; i++)
      {
        cols[0] = pGlobalIdx[k][j][i];
        cols[1] = pGlobalIdx[k][j+1][i];
        countNumNonZeros(cols, 2, lambdaStart, lambdaEnd, d_nnz[localIdx], o_nnz[localIdx]);
        localIdx++;
      }
    }
  }
  // non-zeros for z-momentum equation
  ierr = DMDAGetCorners(wda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  for (k=pstart; k<pstart+p; k++)
  {
    for (j=nstart; j<nstart+n; j++)
    {
      for (i=mstart; i<mstart+m; i++)
      {
        cols[0] = pGlobalIdx[k][j][i];
        cols[1] = pGlobalIdx[k+1][j][i];
        countNumNonZeros(cols, 2, lambdaStart, lambdaEnd, d_nnz[localIdx], o_nnz[localIdx]);
        localIdx++;
      }
    }
  }
  
  // allocate memory for the matrix
  ierr = MatCreate(PETSC_COMM_WORLD, &BNQ); CHKERRQ(ierr);
  ierr = MatSetSizes(BNQ, qLocalSize, lambdaLocalSize, PETSC_DETERMINE, PETSC_DETERMINE); CHKERRQ(ierr);
  ierr = MatSetFromOptions(BNQ); CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(BNQ, 0, d_nnz); CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(BNQ, 0, d_nnz, 0, o_nnz); CHKERRQ(ierr);

  // deallocate d_nnz and o_nnz
  ierr = PetscFree(d_nnz); CHKERRQ(ierr);
  ierr = PetscFree(o_nnz); CHKERRQ(ierr);

  // assemble matrix Q
  localIdx = 0;
  // x-momentum coefficients
  ierr = DMDAGetCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  for (k=pstart; k<pstart+p; k++)
  {
    for (j=nstart; j<nstart+n; j++)
    {
      for (i=mstart; i<mstart+m; i++)
      {
        row = localIdx + qStart;
        cols[0] = pGlobalIdx[k][j][i];
        cols[1] = pGlobalIdx[k][j][i+1];
        ierr = MatSetValues(BNQ, 1, &row, 2, cols, values, INSERT_VALUES); CHKERRQ(ierr);
        localIdx++;
      }
    }
  }
  // y-momentum coefficients
  ierr = DMDAGetCorners(vda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  for (k=pstart; k<pstart+p; k++)
  {
    for (j=nstart; j<nstart+n; j++)
    {
      for (i=mstart; i<mstart+m; i++)
      {
        row = localIdx + qStart;
        cols[0] = pGlobalIdx[k][j][i];
        cols[1] = pGlobalIdx[k][j+1][i];
        ierr = MatSetValues(BNQ, 1, &row, 2, cols, values, INSERT_VALUES); CHKERRQ(ierr);
        localIdx++;
      }
    }
  }
  // z-momentum coefficients
  ierr = DMDAGetCorners(wda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  for (k=pstart; k<pstart+p; k++)
  {
    for (j=nstart; j<nstart+n; j++)
    {
      for (i=mstart; i<mstart+m; i++)
      {
        row = localIdx + qStart;
        cols[0] = pGlobalIdx[k][j][i];
        cols[1] = pGlobalIdx[k+1][j][i];
        ierr = MatSetValues(BNQ, 1, &row, 2, cols, values, INSERT_VALUES); CHKERRQ(ierr);
        localIdx++;
      }
    }
  }
  ierr = DMDAVecRestoreArray(pda, pMapping, &pGlobalIdx); CHKERRQ(ierr);

  ierr = MatAssemblyBegin(BNQ, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(BNQ, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  ierr = MatTranspose(BNQ, MAT_INITIAL_MATRIX, &QT); CHKERRQ(ierr);
  ierr = MatDiagonalScale(BNQ, BN, NULL); CHKERRQ(ierr);

  return 0;
} // generateBNQ