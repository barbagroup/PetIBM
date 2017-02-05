/*! Implementation of the method `generateBNQ` of the class `TairaColoniusSolver`.
 * \file generateBNQ.inl
 */


#include "delta.h"


/*!
 * \brief Assembles the matrices BNQ.
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

  PetscFunctionBeginUser;

  PetscInt i, j, l,        // loop indices
           m, n,           // local number of nodes along each direction
           mstart, nstart; // starting indices
  
  PetscInt localIdx, procIdx;
  PetscInt row, cols[2], BNQ_col, pointIdx;
  PetscReal values[2] = {-1.0, 1.0}; // gradient coefficients
  
  PetscReal value; // to hold the value of the discrete delta function
  PetscReal source[2], // source point, center of the domain of influence
            target[2]; // target point to determine if in domain of influence
  PetscReal disp[2]; // source-target displacement vector
  PetscReal h; // grid-spacing
  // get domain dimensions
  PetscReal widths[2] = {mesh->x[mesh->nx+1] - mesh->x[0],
                         mesh->y[mesh->ny+1] - mesh->y[0]};
  // get boundary types
  BoundaryType bTypes[2] = {flow->boundaries[XPLUS][0].type,
                            flow->boundaries[YPLUS][0].type};
  
  PetscLogEvent  GENERATE_BNQ;
  ierr = PetscLogEventRegister("generateBNQ", 0, &GENERATE_BNQ); CHKERRQ(ierr);
  ierr = PetscLogEventBegin(GENERATE_BNQ, 0, 0, 0, 0); CHKERRQ(ierr);
  
  PetscMPIInt numProcs;
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

  // get ownership range of lambda
  PetscInt lambdaStart, lambdaEnd, lambdaLocalSize;
  ierr = VecGetOwnershipRange(lambda, &lambdaStart, &lambdaEnd); CHKERRQ(ierr);
  lambdaLocalSize = lambdaEnd-lambdaStart;

  // get mapping of pressure values as 2D array
  PetscReal **pMappingArray;
  ierr = DMDAVecGetArray(pda, pMapping, &pMappingArray); CHKERRQ(ierr);

  // determine nnz row by row
  localIdx = 0;
  // rows corresponding to fluxes in x-direction
  ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
  for (j=nstart; j<nstart+n; j++)
  {
    target[1] = 0.5*(mesh->y[j] + mesh->y[j+1]);
    for (i=mstart; i<mstart+m; i++)
    {
      h = mesh->dx[i];
      target[0] = mesh->x[i+1];
      // G portion
      cols[0] = pMappingArray[j][i];
      cols[1] = pMappingArray[j][i+1];
      countNumNonZeros(cols, 2, lambdaStart, lambdaEnd, BNQ_d_nnz[localIdx], BNQ_o_nnz[localIdx]);
      // ET portion
      for (auto &body : bodies)
      { 
        for (l=0; l<body.numPoints; l++)
        {
          pointIdx = body.idxPointsOnProcess[l];
          if (i >= body.I[pointIdx]-2 && i<= body.I[pointIdx]+2 &&
              j >= body.J[pointIdx]-2 && j<= body.J[pointIdx]+2)
          {
            source[0] = body.X[pointIdx];
            source[1] = body.Y[pointIdx];
            if (isInfluenced<2>(target, source, 1.5*h, widths, bTypes, disp))
            {
              BNQ_col = body.globalIdxPoints[l];
              (BNQ_col >= lambdaStart && BNQ_col < lambdaEnd) ? BNQ_d_nnz[localIdx]++ : BNQ_o_nnz[localIdx]++;
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
    target[1] = mesh->y[j+1];
    for (i=mstart; i<mstart+m; i++)
    {
      target[0] = 0.5*(mesh->x[i] + mesh->x[i+1]);
      // G portion
      cols[0] = pMappingArray[j][i];
      cols[1] = pMappingArray[j+1][i];
      countNumNonZeros(cols, 2, lambdaStart, lambdaEnd, BNQ_d_nnz[localIdx], BNQ_o_nnz[localIdx]);
      // ET portion
      for (auto &body : bodies)
      { 
        for (l=0; l<body.numPoints; l++)
        {
          pointIdx = body.idxPointsOnProcess[l];
          if (i >= body.I[pointIdx]-2 && i<= body.I[pointIdx]+2 &&
              j >= body.J[pointIdx]-2 && j<= body.J[pointIdx]+2)
          {
            source[0] = body.X[pointIdx];
            source[1] = body.Y[pointIdx];
            if (isInfluenced<2>(target, source, 1.5*h, widths, bTypes, disp))
            {
              BNQ_col = body.globalIdxPoints[l] + 1;
              (BNQ_col >= lambdaStart && BNQ_col < lambdaEnd) ? BNQ_d_nnz[localIdx]++ : BNQ_o_nnz[localIdx]++;
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

  // deallocate nnz arrays
  // BNQ
  ierr = PetscFree(BNQ_d_nnz); CHKERRQ(ierr);
  ierr = PetscFree(BNQ_o_nnz); CHKERRQ(ierr);

  // assemble matrices Q and ET row by row
  localIdx = 0;
  // rows corresponding to fluxes in x-direction
  ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
  for (j=nstart; j<nstart+n; j++)
  {
    target[1] = 0.5*(mesh->y[j] + mesh->y[j+1]);
    for (i=mstart; i<mstart+m; i++)
    {
      h = mesh->dx[i];
      target[0] = mesh->x[i+1];
      row = localIdx + qStart;
      // G portion
      cols[0] = pMappingArray[j][i];
      cols[1] = pMappingArray[j][i+1];
      ierr = MatSetValues(BNQ, 1, &row, 2, cols, values, INSERT_VALUES); CHKERRQ(ierr);
      // ET portion
      for (auto &body : bodies)
      { 
        for (l=0; l<body.numPoints; l++)
        {
          pointIdx = body.idxPointsOnProcess[l];
          if (i >= body.I[pointIdx]-2 && i<= body.I[pointIdx]+2 &&
              j >= body.J[pointIdx]-2 && j<= body.J[pointIdx]+2)
          {
            source[0] = body.X[pointIdx];
            source[1] = body.Y[pointIdx];
            if (isInfluenced<2>(target, source, 1.5*h, widths, bTypes, disp))
            {
              BNQ_col = body.globalIdxPoints[l];
              value = h*delta(disp[0], disp[1], h);
              ierr = MatSetValue(BNQ, row, BNQ_col, value, INSERT_VALUES); CHKERRQ(ierr);
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
    target[1] = mesh->y[j+1];
    for (i=mstart; i<mstart+m; i++)
    {
      target[0] = 0.5*(mesh->x[i] + mesh->x[i+1]);
      row = localIdx + qStart;
      // G portion
      cols[0] = pMappingArray[j][i];
      cols[1] = pMappingArray[j+1][i];
      ierr = MatSetValues(BNQ, 1, &row, 2, cols, values, INSERT_VALUES); CHKERRQ(ierr);
      // ET portion
      for (auto &body : bodies)
      { 
        for (l=0; l<body.numPoints; l++)
        {
          pointIdx = body.idxPointsOnProcess[l];
          if (i >= body.I[pointIdx]-2 && i<= body.I[pointIdx]+2 &&
              j >= body.J[pointIdx]-2 && j<= body.J[pointIdx]+2)
          {
            source[0] = body.X[pointIdx];
            source[1] = body.Y[pointIdx];
            if (isInfluenced<2>(target, source, 1.5*h, widths, bTypes, disp))
            {
              BNQ_col = body.globalIdxPoints[l] + 1;
              value = h*delta(disp[0], disp[1], h);
              ierr = MatSetValue(BNQ, row, BNQ_col, value, INSERT_VALUES); CHKERRQ(ierr);
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

  ierr = MatTranspose(BNQ, MAT_INITIAL_MATRIX, &QT); CHKERRQ(ierr);
  ierr = MatDiagonalScale(BNQ, BN, NULL); CHKERRQ(ierr);
  
  ierr = PetscLogEventEnd(GENERATE_BNQ, 0, 0, 0, 0); CHKERRQ(ierr);

  ierr = PetscObjectViewFromOptions((PetscObject) BNQ, NULL, "-BNQ_mat_view"); CHKERRQ(ierr);
  ierr = PetscObjectViewFromOptions((PetscObject) QT, NULL, "-QT_mat_view"); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // generateBNQ


// three-dimensional specialization
template <>
PetscErrorCode TairaColoniusSolver<3>::generateBNQ()
{
  PetscErrorCode ierr;

  PetscInt i, j, k, l,             // loop indices
           m, n, p,                // local number of nodes along each direction
           mstart, nstart, pstart; // starting indices

  PetscInt localIdx, procIdx;
  PetscInt row, cols[2], BNQ_col, pointIdx;
  PetscReal values[2] = {-1.0, 1.0}, value;
  
  PetscReal source[3], target[3];
  PetscReal disp[3];
  PetscReal h;
  PetscReal widths[3] = {mesh->x[mesh->nx+1] - mesh->x[0],
                         mesh->y[mesh->ny+1] - mesh->y[0],
                         mesh->z[mesh->nz+1] - mesh->z[0]};
  BoundaryType bTypes[3] = {flow->boundaries[XPLUS][0].type,
                            flow->boundaries[YPLUS][0].type,
                            flow->boundaries[ZPLUS][0].type};
  
  PetscLogEvent  GENERATE_BNQ;
  ierr = PetscLogEventRegister("generateBNQ", 0, &GENERATE_BNQ); CHKERRQ(ierr);
  ierr = PetscLogEventBegin(GENERATE_BNQ, 0, 0, 0, 0); CHKERRQ(ierr);
  
  PetscMPIInt numProcs;
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
  
  // get ownership range of lambda
  PetscInt lambdaStart, lambdaEnd, lambdaLocalSize;
  ierr = VecGetOwnershipRange(lambda, &lambdaStart, &lambdaEnd); CHKERRQ(ierr);
  lambdaLocalSize = lambdaEnd-lambdaStart;

  // get mapping of pressure values
  PetscReal ***pMappingArray;
  ierr = DMDAVecGetArray(pda, pMapping, &pMappingArray); CHKERRQ(ierr);

  // determine nnz row by row
  localIdx = 0;
  // rows corresponding to fluxes in x-direction
  ierr = DMDAGetCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  for (k=pstart; k<pstart+p; k++)
  {
    target[2] = 0.5*(mesh->z[k] + mesh->z[k+1]);
    for (j=nstart; j<nstart+n; j++)
    {
      target[1] = 0.5*(mesh->y[j] + mesh->y[j+1]);
      for (i=mstart; i<mstart+m; i++)
      {
        h = mesh->dx[i];
        target[0] = mesh->x[i+1];
        // G portion
        cols[0] = pMappingArray[k][j][i];
        cols[1] = pMappingArray[k][j][i+1];
        countNumNonZeros(cols, 2, lambdaStart, lambdaEnd, BNQ_d_nnz[localIdx], BNQ_o_nnz[localIdx]);
        // ET portion
        for (auto &body : bodies)
        { 
          for (l=0; l<body.numPoints; l++)
          {
            pointIdx = body.idxPointsOnProcess[l];
            if (i >= body.I[pointIdx]-2 && i<= body.I[pointIdx]+2 &&
                j >= body.J[pointIdx]-2 && j<= body.J[pointIdx]+2 &&
                k >= body.K[pointIdx]-2 && k<= body.K[pointIdx]+2)
            {
              source[0] = body.X[pointIdx];
              source[1] = body.Y[pointIdx];
              source[2] = body.Z[pointIdx];
              if (isInfluenced<3>(target, source, 1.5*h, widths, bTypes, disp))
              {
                BNQ_col = body.globalIdxPoints[l];
                (BNQ_col >= lambdaStart && BNQ_col < lambdaEnd) ? BNQ_d_nnz[localIdx]++ : BNQ_o_nnz[localIdx]++;
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
    target[2] = 0.5*(mesh->z[k] + mesh->z[k+1]);
    for (j=nstart; j<nstart+n; j++)
    {
      h = mesh->dy[j];
      target[1] = mesh->y[j+1];
      for (i=mstart; i<mstart+m; i++)
      {
        target[0] = 0.5*(mesh->x[i] + mesh->x[i+1]);
        // G portion
        cols[0] = pMappingArray[k][j][i];
        cols[1] = pMappingArray[k][j+1][i];
        countNumNonZeros(cols, 2, lambdaStart, lambdaEnd, BNQ_d_nnz[localIdx], BNQ_o_nnz[localIdx]);
        // ET portion
        for (auto &body : bodies)
        { 
          for (l=0; l<body.numPoints; l++)
          {
            pointIdx = body.idxPointsOnProcess[l];
            if (i >= body.I[pointIdx]-2 && i<= body.I[pointIdx]+2 &&
                j >= body.J[pointIdx]-2 && j<= body.J[pointIdx]+2 &&
                k >= body.K[pointIdx]-2 && k<= body.K[pointIdx]+2)
            {
              source[0] = body.X[pointIdx];
              source[1] = body.Y[pointIdx];
              source[2] = body.Z[pointIdx];
              if (isInfluenced<3>(target, source, 1.5*h, widths, bTypes, disp))
              {
                BNQ_col = body.globalIdxPoints[l] + 1;
                (BNQ_col >= lambdaStart && BNQ_col < lambdaEnd) ? BNQ_d_nnz[localIdx]++ : BNQ_o_nnz[localIdx]++;
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
    target[2] = mesh->z[k+1];
    for (j=nstart; j<nstart+n; j++)
    {
      target[1] = 0.5*(mesh->y[j] + mesh->y[j+1]);
      for (i=mstart; i<mstart+m; i++)
      {
        target[0] = 0.5*(mesh->x[i] + mesh->x[i+1]);
        // G portion
        cols[0] = pMappingArray[k][j][i];
        cols[1] = pMappingArray[k+1][j][i];
        countNumNonZeros(cols, 2, lambdaStart, lambdaEnd, BNQ_d_nnz[localIdx], BNQ_o_nnz[localIdx]);
        // ET portion
        for (auto &body : bodies)
        { 
          for (l=0; l<body.numPoints; l++)
          {
            pointIdx = body.idxPointsOnProcess[l];
            if (i >= body.I[pointIdx]-2 && i<= body.I[pointIdx]+2 &&
                j >= body.J[pointIdx]-2 && j<= body.J[pointIdx]+2 &&
                k >= body.K[pointIdx]-2 && k<= body.K[pointIdx]+2)
            {
              source[0] = body.X[pointIdx];
              source[1] = body.Y[pointIdx];
              source[2] = body.Z[pointIdx];
              if (isInfluenced<3>(target, source, 1.5*h, widths, bTypes, disp))
              {
                BNQ_col = body.globalIdxPoints[l] + 2;
                (BNQ_col >= lambdaStart && BNQ_col < lambdaEnd) ? BNQ_d_nnz[localIdx]++ : BNQ_o_nnz[localIdx]++;
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

  // deallocate nnz arrays
  // BNQ
  ierr = PetscFree(BNQ_d_nnz); CHKERRQ(ierr);
  ierr = PetscFree(BNQ_o_nnz); CHKERRQ(ierr);

  // assemble matrix Q row by row
  localIdx = 0;
  // rows corresponding to fluxes in x-direction
  ierr = DMDAGetCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  for (k=pstart; k<pstart+p; k++)
  {
    target[2] = 0.5*(mesh->z[k] + mesh->z[k+1]);
    for (j=nstart; j<nstart+n; j++)
    {
      target[1] = 0.5*(mesh->y[j] + mesh->y[j+1]);
      for (i=mstart; i<mstart+m; i++)
      {
        h = mesh->dx[i];
        target[0] = mesh->x[i+1];
        row = localIdx + qStart;
        // G portion
        cols[0] = pMappingArray[k][j][i];
        cols[1] = pMappingArray[k][j][i+1];
        ierr = MatSetValues(BNQ, 1, &row, 2, cols, values, INSERT_VALUES); CHKERRQ(ierr);
        // ET portion
        for (auto &body : bodies)
        { 
          for (l=0; l<body.numPoints; l++)
          {
            pointIdx = body.idxPointsOnProcess[l];
            if (i >= body.I[pointIdx]-2 && i<= body.I[pointIdx]+2 &&
                j >= body.J[pointIdx]-2 && j<= body.J[pointIdx]+2 &&
                k >= body.K[pointIdx]-2 && k<= body.K[pointIdx]+2)
            {
              source[0] = body.X[pointIdx];
              source[1] = body.Y[pointIdx];
              source[2] = body.Z[pointIdx];
              if (isInfluenced<3>(target, source, 1.5*h, widths, bTypes, disp))
              {
                BNQ_col = body.globalIdxPoints[l];
                value= h*delta(disp[0], disp[1], disp[2], h);
                ierr = MatSetValue(BNQ, row, BNQ_col, value, INSERT_VALUES); CHKERRQ(ierr);
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
    target[2] = 0.5*(mesh->z[k] + mesh->z[k+1]);
    for (j=nstart; j<nstart+n; j++)
    {
      h = mesh->dy[j];
      target[1] = mesh->y[j+1];
      for (i=mstart; i<mstart+m; i++)
      {
        target[0] = 0.5*(mesh->x[i] + mesh->x[i+1]);
        row = localIdx + qStart;
        // G portion
        cols[0] = pMappingArray[k][j][i];
        cols[1] = pMappingArray[k][j+1][i];
        ierr = MatSetValues(BNQ, 1, &row, 2, cols, values, INSERT_VALUES); CHKERRQ(ierr);
        // ET portion
        for (auto &body : bodies)
        { 
          for (l=0; l<body.numPoints; l++)
          {
            pointIdx = body.idxPointsOnProcess[l];
            if (i >= body.I[pointIdx]-2 && i<= body.I[pointIdx]+2 &&
                j >= body.J[pointIdx]-2 && j<= body.J[pointIdx]+2 &&
                k >= body.K[pointIdx]-2 && k<= body.K[pointIdx]+2)
            {
              source[0] = body.X[pointIdx];
              source[1] = body.Y[pointIdx];
              source[2] = body.Z[pointIdx];
              if (isInfluenced<3>(target, source, 1.5*h, widths, bTypes, disp))
              {
                BNQ_col = body.globalIdxPoints[l] + 1;
                value= h*delta(disp[0], disp[1], disp[2], h);
                ierr = MatSetValue(BNQ, row, BNQ_col, value, INSERT_VALUES); CHKERRQ(ierr);
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
    target[2] = mesh->z[k+1];
    for (j=nstart; j<nstart+n; j++)
    {
      target[1] = 0.5*(mesh->y[j] + mesh->y[j+1]);
      for (i=mstart; i<mstart+m; i++)
      {
        target[0] = 0.5*(mesh->x[i] + mesh->x[i+1]);
        row = localIdx + qStart;
        // G portion
        cols[0] = pMappingArray[k][j][i];
        cols[1] = pMappingArray[k+1][j][i];
        ierr = MatSetValues(BNQ, 1, &row, 2, cols, values, INSERT_VALUES); CHKERRQ(ierr);
        // ET portion
        for (auto &body : bodies)
        { 
          for (l=0; l<body.numPoints; l++)
          {
            pointIdx = body.idxPointsOnProcess[l];
            if (i >= body.I[pointIdx]-2 && i<= body.I[pointIdx]+2 &&
                j >= body.J[pointIdx]-2 && j<= body.J[pointIdx]+2 &&
                k >= body.K[pointIdx]-2 && k<= body.K[pointIdx]+2)
            {
              source[0] = body.X[pointIdx];
              source[1] = body.Y[pointIdx];
              source[2] = body.Z[pointIdx];
              if (isInfluenced<3>(target, source, 1.5*h, widths, bTypes, disp))
              {
                BNQ_col = body.globalIdxPoints[l] + 2;
                value= h*delta(disp[0], disp[1], disp[2], h);
                ierr = MatSetValue(BNQ, row, BNQ_col, value, INSERT_VALUES); CHKERRQ(ierr);
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

  ierr = MatTranspose(BNQ, MAT_INITIAL_MATRIX, &QT); CHKERRQ(ierr);
  ierr = MatDiagonalScale(BNQ, BN, NULL); CHKERRQ(ierr);
  
  ierr = PetscLogEventEnd(GENERATE_BNQ, 0, 0, 0, 0); CHKERRQ(ierr);

  ierr = PetscObjectViewFromOptions((PetscObject) BNQ, NULL, "-BNQ_mat_view"); CHKERRQ(ierr);
  ierr = PetscObjectViewFromOptions((PetscObject) QT, NULL, "-QT_mat_view"); CHKERRQ(ierr);

  return 0;
} // generateBNQ
