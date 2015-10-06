/***************************************************************************//**
 * \file generateBNQ.inl
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the method `generateBNQ` of the class `TairaColoniusSolver`.
 */


/**
 * \brief Assembles the matrices BNQ and ET.
 */
template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::generateBNQ()
{
  return 0;
} // generateBNQ


// two-dimensional specialization
template <>
PetscErrorCode TairaColoniusSolver<2>::generateBNQ()
{
  PetscErrorCode ierr;

  PetscInt i, j,           // loop indices
           m, n,           // local number of nodes along each direction
           mstart, nstart; // starting indices
  
  PetscInt localIdx, procIdx;
  PetscInt row, cols[2], BNQ_col, ET_col, value;
  PetscReal values[2] = {-1.0, 1.0};
  PetscReal disp[2];
  PetscReal xCoord, yCoord, h;
  
  PetscLogEvent  GENERATE_BNQ;
  ierr = PetscLogEventRegister("generateBNQ", 0, &GENERATE_BNQ); CHKERRQ(ierr);
  ierr = PetscLogEventBegin(GENERATE_BNQ, 0, 0, 0, 0); CHKERRQ(ierr);
  
  PetscInt numProcs;
  ierr = MPI_Comm_size(PETSC_COMM_WORLD, &numProcs); CHKERRQ(ierr);
  
  // get ownership range of q
  PetscInt qStart, qEnd, qLocalSize;
  ierr = VecGetOwnershipRange(q, &qStart, &qEnd); CHKERRQ(ierr);
  qLocalSize = qEnd-qStart;

  // create arrays to store number of non-zero (nnz) values
  // BNQ
  PetscInt *BNQ_d_nnz, // nnz on diagonal
           *BNQ_o_nnz; // nnz off diagonal
  ierr = PetscMalloc(qLocalSize*sizeof(PetscInt), &BNQ_d_nnz); CHKERRQ(ierr);
  ierr = PetscMalloc(qLocalSize*sizeof(PetscInt), &BNQ_o_nnz); CHKERRQ(ierr);
  // ET
  PetscInt *ET_d_nnz, // nnz on diagonal
           *ET_o_nnz; // nnz off diagonal
  ierr = PetscMalloc(qLocalSize*sizeof(PetscInt), &ET_d_nnz); CHKERRQ(ierr);
  ierr = PetscMalloc(qLocalSize*sizeof(PetscInt), &ET_o_nnz); CHKERRQ(ierr);

  // get ownership range of lambda
  PetscInt lambdaStart, lambdaEnd, lambdaLocalSize;
  ierr = VecGetOwnershipRange(lambda, &lambdaStart, &lambdaEnd); CHKERRQ(ierr);
  lambdaLocalSize = lambdaEnd-lambdaStart;

  // get ownership range of f
  Vec fGlobal;
  ierr = DMCompositeGetAccess(lambdaPack, lambda,  NULL, &fGlobal); CHKERRQ(ierr);
  PetscInt fStart, fEnd, fLocalSize;
  ierr = VecGetOwnershipRange(fGlobal, &fStart, &fEnd); CHKERRQ(ierr);
  fLocalSize = fEnd-fStart;
  ierr = DMCompositeRestoreAccess(lambdaPack, lambda,  NULL, &fGlobal); CHKERRQ(ierr);

  // get mapping of pressure values
  PetscReal **pMappingArray;
  ierr = DMDAVecGetArray(pda, pMapping, &pMappingArray); CHKERRQ(ierr);

  // determine nnz row by row
  localIdx = 0;
  // rows corresponding to fluxes in x-direction
  ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
  for (j=nstart; j<nstart+n; j++)
  {
    yCoord = 0.5*(mesh->y[j] + mesh->y[j+1]);
    for (i=mstart; i<mstart+m; i++)
    {
      h = mesh->dx[i];
      xCoord = mesh->x[i+1];
      // G portion
      cols[0] = pMappingArray[j][i];
      cols[1] = pMappingArray[j][i+1];
      countNumNonZeros(cols, 2, lambdaStart, lambdaEnd, BNQ_d_nnz[localIdx], BNQ_o_nnz[localIdx]);
      // ET portion
      ET_d_nnz[localIdx] = 0;
      ET_o_nnz[localIdx] = 0;
      PetscInt numPhi = 0;
      for (procIdx=0; procIdx<numProcs; procIdx++)
      {
        numPhi += numPhiOnProcess[procIdx];
        for (auto l=boundaryPointIndices[procIdx].begin(); l!=boundaryPointIndices[procIdx].end(); l++)
        {
          if (i >= bodies[0].I[*l]-2 && i <= bodies[0].I[*l]+2 
              && j >= bodies[0].J[*l]-2 && j <= bodies[0].J[*l]+2)
          {
            if (isInfluenced(xCoord, yCoord, bodies[0].X[*l], bodies[0].Y[*l], 1.5*h, disp))
            {
              BNQ_col = globalIndexMapping[*l];
              (BNQ_col >= lambdaStart && BNQ_col < lambdaEnd) ? BNQ_d_nnz[localIdx]++ : BNQ_o_nnz[localIdx]++;
              ET_col = globalIndexMapping[*l] - numPhi;
              (ET_col >= fStart && ET_col < fEnd) ? ET_d_nnz[localIdx]++ : ET_o_nnz[localIdx]++;
            }
          }
        }
      }
      localIdx++;
    }
  }
  // rows corresponding to fluxes in y-direction
  ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
  for (j=nstart; j<nstart+n; j++)
  {
    h = mesh->dy[j];
    yCoord = mesh->y[j+1];
    for (i=mstart; i<mstart+m; i++)
    {
      xCoord = 0.5*(mesh->x[i] + mesh->x[i+1]);
      // G portion
      cols[0] = pMappingArray[j][i];
      cols[1] = pMappingArray[j+1][i];
      countNumNonZeros(cols, 2, lambdaStart, lambdaEnd, BNQ_d_nnz[localIdx], BNQ_o_nnz[localIdx]);
      // ET portion
      ET_d_nnz[localIdx] = 0;
      ET_o_nnz[localIdx] = 0;
      PetscInt numPhi = 0;
      for (procIdx=0; procIdx<numProcs; procIdx++)
      {
        numPhi += numPhiOnProcess[procIdx];
        for (auto l=boundaryPointIndices[procIdx].begin(); l!=boundaryPointIndices[procIdx].end(); l++)
        {
          if (i >= bodies[0].I[*l]-2 && i <= bodies[0].I[*l]+2 
              && j >= bodies[0].J[*l]-2 && j <= bodies[0].J[*l]+2)
          {
            if (isInfluenced(xCoord, yCoord, bodies[0].X[*l], bodies[0].Y[*l], 1.5*h, disp))
            {
              BNQ_col = globalIndexMapping[*l]+1;
              (BNQ_col>=lambdaStart && BNQ_col<lambdaEnd) ? BNQ_d_nnz[localIdx]++ : BNQ_o_nnz[localIdx]++;
              ET_col = globalIndexMapping[*l] - numPhi + 1;
              (ET_col>=fStart && ET_col<fEnd) ? ET_d_nnz[localIdx]++ : ET_o_nnz[localIdx]++;
            }
          }
        }
      }
      localIdx++;
    }
  }
  
  // allocate memory for matrices
  // BNQ
  ierr = MatCreate(PETSC_COMM_WORLD, &BNQ); CHKERRQ(ierr);
  ierr = MatSetSizes(BNQ, qLocalSize, lambdaLocalSize, PETSC_DETERMINE, PETSC_DETERMINE); CHKERRQ(ierr);
  ierr = MatSetFromOptions(BNQ); CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(BNQ, 0, BNQ_d_nnz); CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(BNQ, 0, BNQ_d_nnz, 0, BNQ_o_nnz); CHKERRQ(ierr);
  ierr = MatSetOption(BNQ, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE); CHKERRQ(ierr);
  // ET
  ierr = MatCreate(PETSC_COMM_WORLD, &ET); CHKERRQ(ierr);
  ierr = MatSetSizes(ET, qLocalSize, fLocalSize, PETSC_DETERMINE, PETSC_DETERMINE); CHKERRQ(ierr);
  ierr = MatSetFromOptions(ET); CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(ET, 0, ET_d_nnz); CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(ET, 0, ET_d_nnz, 0, ET_o_nnz); CHKERRQ(ierr);
  ierr = MatSetOption(ET, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE); CHKERRQ(ierr);

  // deallocate nnz arrays
  // BNQ
  ierr = PetscFree(BNQ_d_nnz); CHKERRQ(ierr);
  ierr = PetscFree(BNQ_o_nnz); CHKERRQ(ierr);
  // ET
  ierr = PetscFree(ET_d_nnz); CHKERRQ(ierr);
  ierr = PetscFree(ET_o_nnz); CHKERRQ(ierr);

  // assemble matrices Q and ET row by row
  localIdx = 0;
  // rows corresponding to fluxes in x-direction
  ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
  for (j=nstart; j<nstart+n; j++)
  {
    yCoord = 0.5*(mesh->y[j] + mesh->y[j+1]);
    for (i=mstart; i<mstart+m; i++)
    {
      h = mesh->dx[i];
      xCoord = mesh->x[i+1];
      row = localIdx + qStart;
      // G portion
      cols[0] = pMappingArray[j][i];
      cols[1] = pMappingArray[j][i+1];
      ierr = MatSetValues(BNQ, 1, &row, 2, cols, values, INSERT_VALUES); CHKERRQ(ierr);
      // ET portion
      PetscInt numPhi = 0;
      for (procIdx=0; procIdx<numProcs; procIdx++)
      {
        numPhi += numPhiOnProcess[procIdx];
        for (auto l=boundaryPointIndices[procIdx].begin(); l!=boundaryPointIndices[procIdx].end(); l++)
        {
          if (i >= bodies[0].I[*l]-2 && i <= bodies[0].I[*l]+2 
              && j >= bodies[0].J[*l]-2 && j <= bodies[0].J[*l]+2)
          {
            if (isInfluenced(xCoord, yCoord, bodies[0].X[*l], bodies[0].Y[*l], 1.5*h, disp))
            {
              BNQ_col  = globalIndexMapping[*l];
              value= h*delta(disp[0], disp[1], h);
              ierr = MatSetValue(BNQ, row, BNQ_col, value, INSERT_VALUES); CHKERRQ(ierr);
              ET_col  = globalIndexMapping[*l] - numPhi;
              ierr = MatSetValue(ET, row, ET_col, value, INSERT_VALUES); CHKERRQ(ierr);
            }
          }
        }
      }
      localIdx++;
    }
  }
  // rows corresponding to fluxes in y-direction
  ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
  for (j=nstart; j<nstart+n; j++)
  {
    h = mesh->dy[j];
    yCoord = mesh->y[j+1];
    for (i=mstart; i<mstart+m; i++)
    {
      xCoord = 0.5*(mesh->x[i] + mesh->x[i+1]);
      row = localIdx + qStart;
      // G portion
      cols[0] = pMappingArray[j][i];
      cols[1] = pMappingArray[j+1][i];
      ierr = MatSetValues(BNQ, 1, &row, 2, cols, values, INSERT_VALUES); CHKERRQ(ierr);
      // ET portion
      PetscInt numPhi = 0;
      for (procIdx=0; procIdx<numProcs; procIdx++)
      {
        numPhi += numPhiOnProcess[procIdx];
        for (auto l=boundaryPointIndices[procIdx].begin(); l!=boundaryPointIndices[procIdx].end(); l++)
        {
          if (i >= bodies[0].I[*l]-2 && i <= bodies[0].I[*l]+2 
              && j >= bodies[0].J[*l]-2 && j <= bodies[0].J[*l]+2)
          {
            if (isInfluenced(xCoord, yCoord, bodies[0].X[*l], bodies[0].Y[*l], 1.5*h, disp))
            {
              BNQ_col = globalIndexMapping[*l] + 1;
              value= h*delta(disp[0], disp[1], h);
              ierr = MatSetValue(BNQ, row, BNQ_col, value, INSERT_VALUES); CHKERRQ(ierr);
              ET_col = globalIndexMapping[*l] - numPhi + 1;
              ierr = MatSetValue(ET, row, ET_col, value, INSERT_VALUES); CHKERRQ(ierr);
            }
          }
        }
      }
      localIdx++;
    }
  }
  ierr = DMDAVecRestoreArray(pda, pMapping, &pMappingArray); CHKERRQ(ierr);

  // assemble the matrices
  // BNQ
  ierr = MatAssemblyBegin(BNQ, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(BNQ, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  // ET
  ierr = MatAssemblyBegin(ET, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(ET, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  ierr = MatTranspose(BNQ, MAT_INITIAL_MATRIX, &QT); CHKERRQ(ierr);
  ierr = MatDiagonalScale(BNQ, BN, NULL); CHKERRQ(ierr);
  
  ierr = PetscLogEventEnd(GENERATE_BNQ, 0, 0, 0, 0); CHKERRQ(ierr);

  return 0;
} // generateBNQ


// three-dimensional specialization
template <>
PetscErrorCode TairaColoniusSolver<3>::generateBNQ()
{
  PetscErrorCode ierr;

  PetscInt i, j, k,                // loop indices
           m, n, p,                // local number of nodes along each direction
           mstart, nstart, pstart; // starting indices

  PetscInt localIdx, procIdx;
  PetscInt row, cols[2], BNQ_col, ET_col, value;
  PetscReal values[2] = {-1.0, 1.0};
  PetscReal disp[3];
  PetscReal xCoord, yCoord, zCoord, h;
  
  PetscLogEvent  GENERATE_BNQ;
  ierr = PetscLogEventRegister("generateBNQ", 0, &GENERATE_BNQ); CHKERRQ(ierr);
  ierr = PetscLogEventBegin(GENERATE_BNQ, 0, 0, 0, 0); CHKERRQ(ierr);
  
  PetscInt numProcs;
  ierr = MPI_Comm_size(PETSC_COMM_WORLD, &numProcs); CHKERRQ(ierr);
  
  // get ownership range of q
  PetscInt qStart, qEnd, qLocalSize;
  ierr = VecGetOwnershipRange(q, &qStart, &qEnd); CHKERRQ(ierr);
  qLocalSize = qEnd-qStart;

  // create arrays to store number of non-zero (nnz) values
  // BNQ
  PetscInt *BNQ_d_nnz, // nnz on diagonal
           *BNQ_o_nnz; // nnz off diagonal
  ierr = PetscMalloc(qLocalSize*sizeof(PetscInt), &BNQ_d_nnz); CHKERRQ(ierr);
  ierr = PetscMalloc(qLocalSize*sizeof(PetscInt), &BNQ_o_nnz); CHKERRQ(ierr);
  // ET
  PetscInt *ET_d_nnz, // nnz on diagonal
           *ET_o_nnz; // nnz off diagonal
  ierr = PetscMalloc(qLocalSize*sizeof(PetscInt), &ET_d_nnz); CHKERRQ(ierr);
  ierr = PetscMalloc(qLocalSize*sizeof(PetscInt), &ET_o_nnz); CHKERRQ(ierr);
  
  // get ownership range of lambda
  PetscInt lambdaStart, lambdaEnd, lambdaLocalSize;
  ierr = VecGetOwnershipRange(lambda, &lambdaStart, &lambdaEnd); CHKERRQ(ierr);
  lambdaLocalSize = lambdaEnd-lambdaStart;

  // get ownership range of f
  Vec fGlobal;
  ierr = DMCompositeGetAccess(lambdaPack, lambda, NULL, &fGlobal); CHKERRQ(ierr);
  PetscInt fStart, fEnd, fLocalSize;
  ierr = VecGetOwnershipRange(fGlobal, &fStart, &fEnd); CHKERRQ(ierr);
  fLocalSize = fEnd-fStart;
  ierr = DMCompositeRestoreAccess(lambdaPack, lambda, NULL, &fGlobal); CHKERRQ(ierr);

  // get mapping of pressure values
  PetscReal ***pMappingArray;
  ierr = DMDAVecGetArray(pda, pMapping, &pMappingArray); CHKERRQ(ierr);

  // determine nnz row by row
  localIdx = 0;
  // rows corresponding to fluxes in x-direction
  ierr = DMDAGetCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  for (k=pstart; k<pstart+p; k++)
  {
    zCoord = 0.5*(mesh->z[k] + mesh->z[k+1]);
    for (j=nstart; j<nstart+n; j++)
    {
      yCoord = 0.5*(mesh->y[j] + mesh->y[j+1]);
      for (i=mstart; i<mstart+m; i++)
      {
        h = mesh->dx[i];
        xCoord = mesh->x[i+1];
        // G portion
        cols[0] = pMappingArray[k][j][i];
        cols[1] = pMappingArray[k][j][i+1];
        countNumNonZeros(cols, 2, lambdaStart, lambdaEnd, BNQ_d_nnz[localIdx], BNQ_o_nnz[localIdx]);
        // ET portion
        ET_d_nnz[localIdx] = 0;
        ET_o_nnz[localIdx] = 0;
        PetscInt numPhi = 0;
        for (procIdx=0; procIdx<numProcs; procIdx++)
        {
          numPhi += numPhiOnProcess[procIdx];
          for (auto l=boundaryPointIndices[procIdx].begin(); l!=boundaryPointIndices[procIdx].end(); l++)
          {
            if (i >= bodies[0].I[*l]-2 && i <= bodies[0].I[*l]+2 
                && j >= bodies[0].J[*l]-2 && j <= bodies[0].J[*l]+2 
                && k >= bodies[0].K[*l]-2 && k <= bodies[0].K[*l]+2)
            {
              if (isInfluenced(xCoord, yCoord, zCoord, bodies[0].X[*l], bodies[0].Y[*l], bodies[0].Z[*l], 1.5*h, disp))
              {
                BNQ_col = globalIndexMapping[*l];
                (BNQ_col >= lambdaStart && BNQ_col < lambdaEnd) ? BNQ_d_nnz[localIdx]++ : BNQ_o_nnz[localIdx]++;
                ET_col = globalIndexMapping[*l] - numPhi;
                (ET_col >= fStart && ET_col < fEnd) ? ET_d_nnz[localIdx]++ : ET_o_nnz[localIdx]++;
              }
            }
          }
        }
        localIdx++;
      }
    }
  }
  // rows corresponding to fluxes in y-direction
  ierr = DMDAGetCorners(vda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  for (k=pstart; k<pstart+p; k++)
  {
    zCoord = 0.5*(mesh->z[k] + mesh->z[k+1]);
    for (j=nstart; j<nstart+n; j++)
    {
      h = mesh->dy[j];
      yCoord = mesh->y[j+1];
      for (i=mstart; i<mstart+m; i++)
      {
        xCoord = 0.5*(mesh->x[i] + mesh->x[i+1]);
        // G portion
        cols[0] = pMappingArray[k][j][i];
        cols[1] = pMappingArray[k][j+1][i];
        countNumNonZeros(cols, 2, lambdaStart, lambdaEnd, BNQ_d_nnz[localIdx], BNQ_o_nnz[localIdx]);
        // ET portion
        ET_d_nnz[localIdx] = 0;
        ET_o_nnz[localIdx] = 0;
        PetscInt numPhi = 0;
        for (procIdx=0; procIdx<numProcs; procIdx++)
        {
          numPhi += numPhiOnProcess[procIdx];
          for (auto l=boundaryPointIndices[procIdx].begin(); l!=boundaryPointIndices[procIdx].end(); l++)
          {
            if (i >= bodies[0].I[*l]-2 && i <= bodies[0].I[*l]+2 
                && j >= bodies[0].J[*l]-2 && j <= bodies[0].J[*l]+2 
                && k >= bodies[0].K[*l]-2 && k <= bodies[0].K[*l]+2)
            {
              if (isInfluenced(xCoord, yCoord, zCoord, bodies[0].X[*l], bodies[0].Y[*l], bodies[0].Z[*l], 1.5*h, disp))
              {
                BNQ_col = globalIndexMapping[*l] + 1;
                (BNQ_col >= lambdaStart && BNQ_col < lambdaEnd) ? BNQ_d_nnz[localIdx]++ : BNQ_o_nnz[localIdx]++;
                ET_col = globalIndexMapping[*l] - numPhi + 1;
                (ET_col >= fStart && ET_col < fEnd) ? ET_d_nnz[localIdx]++ : ET_o_nnz[localIdx]++;
              }
            }
          }
        }
        localIdx++;
      }
    }
  }
  // rows corresponding to fluxes in z-direction
  ierr = DMDAGetCorners(wda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  for (k=pstart; k<pstart+p; k++)
  {
    h = mesh->dz[k];
    zCoord = mesh->z[k+1];
    for (j=nstart; j<nstart+n; j++)
    {
      yCoord = 0.5*(mesh->y[j] + mesh->y[j+1]);
      for (i=mstart; i<mstart+m; i++)
      {
        xCoord = 0.5*(mesh->x[i] + mesh->x[i+1]);
        // G portion
        cols[0] = pMappingArray[k][j][i];
        cols[1] = pMappingArray[k+1][j][i];
        countNumNonZeros(cols, 2, lambdaStart, lambdaEnd, BNQ_d_nnz[localIdx], BNQ_o_nnz[localIdx]);
        // ET portion
        ET_d_nnz[localIdx] = 0;
        ET_o_nnz[localIdx] = 0;
        PetscInt numPhi = 0;
        for (procIdx=0; procIdx<numProcs; procIdx++)
        {
          numPhi += numPhiOnProcess[procIdx];
          for (auto l=boundaryPointIndices[procIdx].begin(); l!=boundaryPointIndices[procIdx].end(); l++)
          {
            if (i >= bodies[0].I[*l]-2 && i <= bodies[0].I[*l]+2 
                && j >= bodies[0].J[*l]-2 && j <= bodies[0].J[*l]+2 
                && k >= bodies[0].K[*l]-2 && k <= bodies[0].K[*l]+2)
            {
              if (isInfluenced(xCoord, yCoord, zCoord, bodies[0].X[*l], bodies[0].Y[*l], bodies[0].Z[*l], 1.5*h, disp))
              {
                BNQ_col = globalIndexMapping[*l] + 2;
                (BNQ_col >= lambdaStart && BNQ_col < lambdaEnd) ? BNQ_d_nnz[localIdx]++ : BNQ_o_nnz[localIdx]++;
                ET_col = globalIndexMapping[*l] - numPhi + 2;
                (ET_col >= fStart && ET_col < fEnd) ? ET_d_nnz[localIdx]++ : ET_o_nnz[localIdx]++;
              }
            }
          }
        }
        localIdx++;
      }
    }
  }
  
  // allocate memory for matrices
  // BNQ
  ierr = MatCreate(PETSC_COMM_WORLD, &BNQ); CHKERRQ(ierr);
  ierr = MatSetSizes(BNQ, qLocalSize, lambdaLocalSize, PETSC_DETERMINE, PETSC_DETERMINE); CHKERRQ(ierr);
  ierr = MatSetFromOptions(BNQ); CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(BNQ, 0, BNQ_d_nnz); CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(BNQ, 0, BNQ_d_nnz, 0, BNQ_o_nnz); CHKERRQ(ierr);
  ierr = MatSetOption(BNQ, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE); CHKERRQ(ierr);
  // ET
  ierr = MatCreate(PETSC_COMM_WORLD, &ET); CHKERRQ(ierr);
  ierr = MatSetSizes(ET, qLocalSize, fLocalSize, PETSC_DETERMINE, PETSC_DETERMINE); CHKERRQ(ierr);
  ierr = MatSetFromOptions(ET); CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(ET, 0, ET_d_nnz); CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(ET, 0, ET_d_nnz, 0, ET_o_nnz); CHKERRQ(ierr);
  ierr = MatSetOption(ET, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE); CHKERRQ(ierr);

  // deallocate nnz arrays
  // BNQ
  ierr = PetscFree(BNQ_d_nnz); CHKERRQ(ierr);
  ierr = PetscFree(BNQ_o_nnz); CHKERRQ(ierr);
  // ET
  ierr = PetscFree(ET_d_nnz); CHKERRQ(ierr);
  ierr = PetscFree(ET_o_nnz); CHKERRQ(ierr);

  // assemble matrix Q row by row
  localIdx = 0;
  // rows corresponding to fluxes in x-direction
  ierr = DMDAGetCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  for (k=pstart; k<pstart+p; k++)
  {
    zCoord = 0.5*(mesh->z[k] + mesh->z[k+1]);
    for (j=nstart; j<nstart+n; j++)
    {
      yCoord = 0.5*(mesh->y[j] + mesh->y[j+1]);
      for (i=mstart; i<mstart+m; i++)
      {
        h = mesh->dx[i];
        xCoord = mesh->x[i+1];
        row = localIdx + qStart;
        // G portion
        cols[0] = pMappingArray[k][j][i];
        cols[1] = pMappingArray[k][j][i+1];
        ierr = MatSetValues(BNQ, 1, &row, 2, cols, values, INSERT_VALUES); CHKERRQ(ierr);
        // ET portion
        PetscInt numPhi = 0;
        for (procIdx=0; procIdx<numProcs; procIdx++)
        {
          numPhi += numPhiOnProcess[procIdx];
          for (auto l=boundaryPointIndices[procIdx].begin(); l!=boundaryPointIndices[procIdx].end(); l++)
          {
            if (i >= bodies[0].I[*l]-2 && i <= bodies[0].I[*l]+2 
                && j >= bodies[0].J[*l]-2 && j <= bodies[0].J[*l]+2 
                && k >= bodies[0].K[*l]-2 && k <= bodies[0].K[*l]+2)
            {
              if (isInfluenced(xCoord, yCoord, zCoord, bodies[0].X[*l], bodies[0].Y[*l], bodies[0].Z[*l], 1.5*h, disp))
              {
                BNQ_col = globalIndexMapping[*l];
                value= h*delta(disp[0], disp[1], disp[2], h);
                ierr = MatSetValue(BNQ, row, BNQ_col, value, INSERT_VALUES); CHKERRQ(ierr);
                ET_col = globalIndexMapping[*l] - numPhi;
                ierr = MatSetValue(ET, row, ET_col, value, INSERT_VALUES); CHKERRQ(ierr);
              }
            }
          }
        }
        localIdx++;
      }
    }
  }
  // rows corresponding to fluxes in y-direction
  ierr = DMDAGetCorners(vda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  for (k=pstart; k<pstart+p; k++)
  {
    zCoord = 0.5*(mesh->z[k] + mesh->z[k+1]);
    for (j=nstart; j<nstart+n; j++)
    {
      h = mesh->dy[j];
      yCoord = mesh->y[j+1];
      for (i=mstart; i<mstart+m; i++)
      {
        xCoord = 0.5*(mesh->x[i] + mesh->x[i+1]);
        row = localIdx + qStart;
        // G portion
        cols[0] = pMappingArray[k][j][i];
        cols[1] = pMappingArray[k][j+1][i];
        ierr = MatSetValues(BNQ, 1, &row, 2, cols, values, INSERT_VALUES); CHKERRQ(ierr);
        // ET portion
        PetscInt numPhi = 0;
        for (procIdx=0; procIdx<numProcs; procIdx++)
        {
          numPhi += numPhiOnProcess[procIdx];
          for (auto l=boundaryPointIndices[procIdx].begin(); l!=boundaryPointIndices[procIdx].end(); l++)
          {
            if (i >= bodies[0].I[*l]-2 && i <= bodies[0].I[*l]+2 
                && j >= bodies[0].J[*l]-2 && j <= bodies[0].J[*l]+2 
                && k >= bodies[0].K[*l]-2 && k <= bodies[0].K[*l]+2)
            {
              if (isInfluenced(xCoord, yCoord, zCoord, bodies[0].X[*l], bodies[0].Y[*l], bodies[0].Z[*l], 1.5*h, disp))
              {
                BNQ_col = globalIndexMapping[*l] + 1;
                value = h*delta(disp[0], disp[1], disp[2], h);
                ierr = MatSetValue(BNQ, row, BNQ_col, value, INSERT_VALUES); CHKERRQ(ierr);
                ET_col = globalIndexMapping[*l] - numPhi + 1;
                ierr = MatSetValue(ET, row, ET_col, value, INSERT_VALUES); CHKERRQ(ierr);
              }
            }
          }
        }
        localIdx++;
      }
    }
  }
  // rows corresponding to fluxes in z-direction
  ierr = DMDAGetCorners(wda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  for (k=pstart; k<pstart+p; k++)
  {
    h = mesh->dz[k];
    zCoord = mesh->z[k+1];
    for (j=nstart; j<nstart+n; j++)
    {
      yCoord = 0.5*(mesh->y[j] + mesh->y[j+1]);
      for (i=mstart; i<mstart+m; i++)
      {
        xCoord = 0.5*(mesh->x[i] + mesh->x[i+1]);
        row = localIdx + qStart;
        // G portion
        cols[0] = pMappingArray[k][j][i];
        cols[1] = pMappingArray[k+1][j][i];
        ierr = MatSetValues(BNQ, 1, &row, 2, cols, values, INSERT_VALUES); CHKERRQ(ierr);
        // ET portion
        PetscInt numPhi = 0;
        for (procIdx=0; procIdx<numProcs; procIdx++)
        {
          numPhi += numPhiOnProcess[procIdx];
          for (auto l=boundaryPointIndices[procIdx].begin(); l!=boundaryPointIndices[procIdx].end(); l++)
          {
            if (i >= bodies[0].I[*l]-2 && i <= bodies[0].I[*l]+2 
                && j >= bodies[0].J[*l]-2 && j <= bodies[0].J[*l]+2 
                && k >= bodies[0].K[*l]-2 && k <= bodies[0].K[*l]+2)
            {
              if(isInfluenced(xCoord, yCoord, zCoord, bodies[0].X[*l], bodies[0].Y[*l], bodies[0].Z[*l], 1.5*h, disp))
              {
                BNQ_col = globalIndexMapping[*l] + 2;
                value = h*delta(disp[0], disp[1], disp[2], h);
                ierr = MatSetValue(BNQ, row, BNQ_col, value, INSERT_VALUES); CHKERRQ(ierr);
                ET_col = globalIndexMapping[*l] - numPhi + 2;
                ierr = MatSetValue(ET, row, ET_col, value, INSERT_VALUES); CHKERRQ(ierr);
              }
            }
          }
        }
        localIdx++;
      }
    }
  }
  ierr = DMDAVecRestoreArray(pda, pMapping, &pMappingArray); CHKERRQ(ierr);

  // assembles matrices
  // BNQ
  ierr = MatAssemblyBegin(BNQ, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(BNQ, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  // ET
  ierr = MatAssemblyBegin(ET, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(ET, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  ierr = MatTranspose(BNQ, MAT_INITIAL_MATRIX, &QT); CHKERRQ(ierr);
  ierr = MatDiagonalScale(BNQ, BN, NULL); CHKERRQ(ierr);
  
  ierr = PetscLogEventEnd(GENERATE_BNQ, 0, 0, 0, 0); CHKERRQ(ierr);

  return 0;
} // generateBNQ
