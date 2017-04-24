/*! Implementation of the method `generateET` of the class `LiEtAlSolver`.
 * \file generateET.inl
 */


#include "delta.h"


/*!
 * \brief Assembles the matrices ET.
 */
template <PetscInt dim>
PetscErrorCode LiEtAlSolver<dim>::generateET()
{
  return 0;
} // generateET


// two-dimensional specialization
template <>
PetscErrorCode LiEtAlSolver<2>::generateET()
{
  PetscErrorCode ierr;

  PetscInt i, j, l,        // loop indices
           m, n,           // local number of nodes along each direction
           mstart, nstart; // starting indices
  
  PetscInt localIdx, procIdx;
  PetscInt row, ET_col, pointIdx;
  
  PetscReal value; // to hold the value of the discrete delta function
  PetscReal source[2], // source point, center of the domain of influence
            target[2]; // target point to determine if in domain of influence
  PetscReal disp[2]; // source-target displacement vector
  PetscReal hx, hy; // grid-spacings
  PetscReal maxDisp[2]; // lengths of the support of the discrete delta function
  // get domain dimensions
  PetscReal widths[2] = {mesh->x[mesh->nx] - mesh->x[0],
                         mesh->y[mesh->ny] - mesh->y[0]};
  // get boundary types
  BoundaryType bTypes[2] = {flow->boundaries[XPLUS][0].type,
                            flow->boundaries[YPLUS][0].type};

  PetscFunctionBeginUser;

  // get ownership range of q
  PetscInt qStart, qEnd, qLocalSize;
  ierr = VecGetOwnershipRange(q, &qStart, &qEnd); CHKERRQ(ierr);
  qLocalSize = qEnd-qStart;

  PetscInt *ET_d_nnz, // nnz on diagonal
           *ET_o_nnz; // nnz off diagonal
  ierr = PetscMalloc(qLocalSize*sizeof(PetscInt), &ET_d_nnz); CHKERRQ(ierr);
  ierr = PetscMalloc(qLocalSize*sizeof(PetscInt), &ET_o_nnz); CHKERRQ(ierr);

  // get ownership range of fTilde
  PetscInt fStart, fEnd, fLocalSize;
  ierr = VecGetOwnershipRange(fTilde, &fStart, &fEnd); CHKERRQ(ierr);
  fLocalSize = fEnd-fStart;

  // determine nnz row by row
  localIdx = 0;
  // rows corresponding to fluxes in x-direction
  ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
  for (j=nstart; j<nstart+n; j++)
  {
    hy = mesh->dy[j];
    target[1] = 0.5*(mesh->y[j] + mesh->y[j+1]);
    maxDisp[1] = 1.5*hy;
    for (i=mstart; i<mstart+m; i++)
    {
      ET_d_nnz[localIdx] = 0;
      ET_o_nnz[localIdx] = 0;
      hx = mesh->dx[i];
      target[0] = mesh->x[i+1];
      maxDisp[0] = 1.5*hx;
      // ET portion
      for (auto &body : bodies)
      { 
        for (l=0; l<body.numPoints; l++)
        {
          pointIdx = body.idxPointsOnProcess[l];
          source[0] = body.X[pointIdx];
          source[1] = body.Y[pointIdx];
          if (isInfluenced<2>(target, source, maxDisp, widths, bTypes, disp))
          {
            ET_col = body.globalIdxPoints[l];
            (ET_col >= fStart && ET_col < fEnd) ? ET_d_nnz[localIdx]++ : ET_o_nnz[localIdx]++;
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
    hy = mesh->dy[j];
    target[1] = mesh->y[j+1];
    maxDisp[1] = 1.5*hy;
    for (i=mstart; i<mstart+m; i++)
    {
      ET_d_nnz[localIdx] = 0;
      ET_o_nnz[localIdx] = 0;
      hx = mesh->dx[i];
      target[0] = 0.5*(mesh->x[i] + mesh->x[i+1]);
      maxDisp[0] = 1.5*hx;
      // ET portion
      for (auto &body : bodies)
      { 
        for (l=0; l<body.numPoints; l++)
        {
          pointIdx = body.idxPointsOnProcess[l];
          source[0] = body.X[pointIdx];
          source[1] = body.Y[pointIdx];
          if (isInfluenced<2>(target, source, maxDisp, widths, bTypes, disp))
          {
            ET_col = body.globalIdxPoints[l] + 1;
            (ET_col >= fStart && ET_col < fEnd) ? ET_d_nnz[localIdx]++ : ET_o_nnz[localIdx]++;
          }
        }
      }
      localIdx++;
    }
  }

  // allocate memory for matrix ET
  ierr = MatCreate(PETSC_COMM_WORLD, &ET); CHKERRQ(ierr);
  ierr = MatSetSizes(ET, qLocalSize, fLocalSize, PETSC_DETERMINE, PETSC_DETERMINE); CHKERRQ(ierr);
  ierr = MatSetFromOptions(ET); CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(ET, 0, ET_d_nnz); CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(ET, 0, ET_d_nnz, 0, ET_o_nnz); CHKERRQ(ierr);
  ierr = MatSetOption(ET, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE); CHKERRQ(ierr);

  // deallocate nnz arrays
  ierr = PetscFree(ET_d_nnz); CHKERRQ(ierr);
  ierr = PetscFree(ET_o_nnz); CHKERRQ(ierr);

  // assemble matrix ET row by row
  localIdx = 0;
  // rows corresponding to fluxes in x-direction
  ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
  for (j=nstart; j<nstart+n; j++)
  {
    hy = mesh->dy[j];
    target[1] = 0.5*(mesh->y[j] + mesh->y[j+1]);
    maxDisp[1] = 1.5*hy;
    for (i=mstart; i<mstart+m; i++)
    {
      hx = mesh->dx[i];
      target[0] = mesh->x[i+1];
      maxDisp[0] = 1.5*hx;
      row = localIdx + qStart;
      // ET portion
      for (auto &body : bodies)
      { 
        for (l=0; l<body.numPoints; l++)
        {
          pointIdx = body.idxPointsOnProcess[l];
          source[0] = body.X[pointIdx];
          source[1] = body.Y[pointIdx];
          if (isInfluenced<2>(target, source, maxDisp, widths, bTypes, disp))
          {
            ET_col = body.globalIdxPoints[l];
            value = hx*delta(disp[0], disp[1], hx, hy);
            ierr = MatSetValue(ET, row, ET_col, value, INSERT_VALUES); CHKERRQ(ierr);
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
    hy = mesh->dy[j];
    target[1] = mesh->y[j+1];
    maxDisp[1] = 1.5*hy;
    for (i=mstart; i<mstart+m; i++)
    {
      hx = mesh->dx[i];
      target[0] = 0.5*(mesh->x[i] + mesh->x[i+1]);
      maxDisp[0] = 1.5*hx;
      row = localIdx + qStart;
      // ET portion
      for (auto &body : bodies)
      { 
        for (l=0; l<body.numPoints; l++)
        {
          pointIdx = body.idxPointsOnProcess[l];
          source[0] = body.X[pointIdx];
          source[1] = body.Y[pointIdx];
          if (isInfluenced<2>(target, source, maxDisp, widths, bTypes, disp))
          {
            ET_col = body.globalIdxPoints[l] + 1;
            value = hy*delta(disp[0], disp[1], hx, hy);
            ierr = MatSetValue(ET, row, ET_col, value, INSERT_VALUES); CHKERRQ(ierr);
          }
        }
      }
      localIdx++;
    }
  }

  ierr = MatAssemblyBegin(ET, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(ET, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  ierr = PetscObjectViewFromOptions((PetscObject) ET, NULL, "-ET_mat_view"); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // generateET


// three-dimensional specialization
template <>
PetscErrorCode LiEtAlSolver<3>::generateET()
{
  PetscErrorCode ierr;

  PetscInt i, j, k, l,             // loop indices
           m, n, p,                // local number of nodes along each direction
           mstart, nstart, pstart; // starting indices

  PetscInt localIdx, procIdx;
  PetscInt row, ET_col, pointIdx;
  PetscReal value;
  
  PetscReal source[3], target[3];
  PetscReal disp[3];
  PetscReal hx, hy, hz; // grid-spacings
  PetscReal maxDisp[3]; // lengths of the support of the discrete delta function
  PetscReal widths[3] = {mesh->x[mesh->nx] - mesh->x[0],
                         mesh->y[mesh->ny] - mesh->y[0],
                         mesh->z[mesh->nz] - mesh->z[0]};
  BoundaryType bTypes[3] = {flow->boundaries[XPLUS][0].type,
                            flow->boundaries[YPLUS][0].type,
                            flow->boundaries[ZPLUS][0].type};

  PetscFunctionBeginUser;

  // get ownership range of q
  PetscInt qStart, qEnd, qLocalSize;
  ierr = VecGetOwnershipRange(q, &qStart, &qEnd); CHKERRQ(ierr);
  qLocalSize = qEnd-qStart;

  // create arrays to store number of non-zero (nnz) values
  PetscInt *ET_d_nnz, // nnz on diagonal
           *ET_o_nnz; // nnz off diagonal
  ierr = PetscMalloc(qLocalSize*sizeof(PetscInt), &ET_d_nnz); CHKERRQ(ierr);
  ierr = PetscMalloc(qLocalSize*sizeof(PetscInt), &ET_o_nnz); CHKERRQ(ierr);

  // get ownership range of fTilde
  PetscInt fStart, fEnd, fLocalSize;
  ierr = VecGetOwnershipRange(fTilde, &fStart, &fEnd); CHKERRQ(ierr);
  fLocalSize = fEnd-fStart;

  // determine nnz row by row
  localIdx = 0;
  // rows corresponding to fluxes in x-direction
  ierr = DMDAGetCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  for (k=pstart; k<pstart+p; k++)
  {
    hz = mesh->dz[k];
    target[2] = 0.5*(mesh->z[k] + mesh->z[k+1]);
    maxDisp[2] = 1.5*hz;
    for (j=nstart; j<nstart+n; j++)
    {
      hy = mesh->dy[j];
      target[1] = 0.5*(mesh->y[j] + mesh->y[j+1]);
      maxDisp[1] = 1.5*hy;
      for (i=mstart; i<mstart+m; i++)
      {
        ET_d_nnz[localIdx] = 0;
        ET_o_nnz[localIdx] = 0;
        hx = mesh->dx[i];
        target[0] = mesh->x[i+1];
        maxDisp[0] = 1.5*hx;
        // ET portion
        for (auto &body : bodies)
        { 
          for (l=0; l<body.numPoints; l++)
          {
            pointIdx = body.idxPointsOnProcess[l];
            source[0] = body.X[pointIdx];
            source[1] = body.Y[pointIdx];
            source[2] = body.Z[pointIdx];
            if (isInfluenced<3>(target, source, maxDisp, widths, bTypes, disp))
            {
              ET_col = body.globalIdxPoints[l];
              (ET_col >= fStart && ET_col < fEnd) ? ET_d_nnz[localIdx]++ : ET_o_nnz[localIdx]++;
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
    hz = mesh->dz[k];
    target[2] = 0.5*(mesh->z[k] + mesh->z[k+1]);
    maxDisp[2] = 1.5*hz;
    for (j=nstart; j<nstart+n; j++)
    {
      hy = mesh->dy[j];
      target[1] = mesh->y[j+1];
      maxDisp[1] = 1.5*hy;
      for (i=mstart; i<mstart+m; i++)
      {
        ET_d_nnz[localIdx] = 0;
        ET_o_nnz[localIdx] = 0;
        hx = mesh->dx[i];
        target[0] = 0.5*(mesh->x[i] + mesh->x[i+1]);
        maxDisp[0] = 1.5*hx;
        // ET portion
        for (auto &body : bodies)
        { 
          for (l=0; l<body.numPoints; l++)
          {
            pointIdx = body.idxPointsOnProcess[l];
            source[0] = body.X[pointIdx];
            source[1] = body.Y[pointIdx];
            source[2] = body.Z[pointIdx];
            if (isInfluenced<3>(target, source, maxDisp, widths, bTypes, disp))
            {
              ET_col = body.globalIdxPoints[l] + 1;
              (ET_col >= fStart && ET_col < fEnd) ? ET_d_nnz[localIdx]++ : ET_o_nnz[localIdx]++;
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
    hz = mesh->dz[k];
    target[2] = mesh->z[k+1];
    maxDisp[2] = 1.5*hz;
    for (j=nstart; j<nstart+n; j++)
    {
      hy = mesh->dy[j];
      target[1] = 0.5*(mesh->y[j] + mesh->y[j+1]);
      maxDisp[1] = 1.5*hy;
      for (i=mstart; i<mstart+m; i++)
      {
        ET_d_nnz[localIdx] = 0;
        ET_o_nnz[localIdx] = 0;
        hx = mesh->dx[i];
        target[0] = 0.5*(mesh->x[i] + mesh->x[i+1]);
        maxDisp[0] = 1.5*hx;
        // ET portion
        for (auto &body : bodies)
        { 
          for (l=0; l<body.numPoints; l++)
          {
            pointIdx = body.idxPointsOnProcess[l];
            source[0] = body.X[pointIdx];
            source[1] = body.Y[pointIdx];
            source[2] = body.Z[pointIdx];
            if (isInfluenced<3>(target, source, maxDisp, widths, bTypes, disp))
            {
              ET_col = body.globalIdxPoints[l] + 2;
              (ET_col >= fStart && ET_col < fEnd) ? ET_d_nnz[localIdx]++ : ET_o_nnz[localIdx]++;
            }
          }
        }
        localIdx++;
      }
    }
  }
  
  // allocate memory for matrix ET
  ierr = MatCreate(PETSC_COMM_WORLD, &ET); CHKERRQ(ierr);
  ierr = MatSetSizes(ET, qLocalSize, fLocalSize, PETSC_DETERMINE, PETSC_DETERMINE); CHKERRQ(ierr);
  ierr = MatSetFromOptions(ET); CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(ET, 0, ET_d_nnz); CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(ET, 0, ET_d_nnz, 0, ET_o_nnz); CHKERRQ(ierr);
  ierr = MatSetOption(ET, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE); CHKERRQ(ierr);

  // deallocate nnz arrays
  ierr = PetscFree(ET_d_nnz); CHKERRQ(ierr);
  ierr = PetscFree(ET_o_nnz); CHKERRQ(ierr);

  // assemble matrix ET row by row
  localIdx = 0;
  // rows corresponding to fluxes in x-direction
  ierr = DMDAGetCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  for (k=pstart; k<pstart+p; k++)
  {
    hz = mesh->dz[k];
    target[2] = 0.5*(mesh->z[k] + mesh->z[k+1]);
    maxDisp[2] = 1.5*hz;
    for (j=nstart; j<nstart+n; j++)
    {
      hy = mesh->dy[j];
      target[1] = 0.5*(mesh->y[j] + mesh->y[j+1]);
      maxDisp[1] = 1.5*hy;
      for (i=mstart; i<mstart+m; i++)
      {
        hx = mesh->dx[i];
        target[0] = mesh->x[i+1];
        maxDisp[0] = 1.5*hx;
        row = localIdx + qStart;
        // ET portion
        for (auto &body : bodies)
        { 
          for (l=0; l<body.numPoints; l++)
          {
            pointIdx = body.idxPointsOnProcess[l];
            source[0] = body.X[pointIdx];
            source[1] = body.Y[pointIdx];
            source[2] = body.Z[pointIdx];
            if (isInfluenced<3>(target, source, maxDisp, widths, bTypes, disp))
            {
              ET_col = body.globalIdxPoints[l];
              value = hx*delta(disp[0], disp[1], disp[2], hx, hy, hz);
              ierr = MatSetValue(ET, row, ET_col, value, INSERT_VALUES); CHKERRQ(ierr);
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
    hz = mesh->dz[k];
    target[2] = 0.5*(mesh->z[k] + mesh->z[k+1]);
    maxDisp[2] = 1.5*hz;
    for (j=nstart; j<nstart+n; j++)
    {
      hy = mesh->dy[j];
      target[1] = mesh->y[j+1];
      maxDisp[1] = 1.5*hy;
      for (i=mstart; i<mstart+m; i++)
      {
        hx = mesh->dx[i];
        target[0] = 0.5*(mesh->x[i] + mesh->x[i+1]);
        maxDisp[0] = 1.5*hx;
        row = localIdx + qStart;
        // ET portion
        for (auto &body : bodies)
        { 
          for (l=0; l<body.numPoints; l++)
          {
            pointIdx = body.idxPointsOnProcess[l];
            source[0] = body.X[pointIdx];
            source[1] = body.Y[pointIdx];
            source[2] = body.Z[pointIdx];
            if (isInfluenced<3>(target, source, maxDisp, widths, bTypes, disp))
            {
              ET_col = body.globalIdxPoints[l] + 1;
              value = hy*delta(disp[0], disp[1], disp[2], hx, hy, hz);
              ierr = MatSetValue(ET, row, ET_col, value, INSERT_VALUES); CHKERRQ(ierr);
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
    hz = mesh->dz[k];
    target[2] = mesh->z[k+1];
    maxDisp[2] = 1.5*hz;
    for (j=nstart; j<nstart+n; j++)
    {
      hy = mesh->dy[j];
      target[1] = 0.5*(mesh->y[j] + mesh->y[j+1]);
      maxDisp[1] = 1.5*hy;
      for (i=mstart; i<mstart+m; i++)
      {
        hx = mesh->dx[i];
        target[0] = 0.5*(mesh->x[i] + mesh->x[i+1]);
        maxDisp[0] = 1.5*hx;
        row = localIdx + qStart;
        // ET portion
        for (auto &body : bodies)
        { 
          for (l=0; l<body.numPoints; l++)
          {
            pointIdx = body.idxPointsOnProcess[l];
            source[0] = body.X[pointIdx];
            source[1] = body.Y[pointIdx];
            source[2] = body.Z[pointIdx];
            if (isInfluenced<3>(target, source, maxDisp, widths, bTypes, disp))
            {
              ET_col = body.globalIdxPoints[l] + 2;
              value = hz*delta(disp[0], disp[1], disp[2], hx, hy, hz);
              ierr = MatSetValue(ET, row, ET_col, value, INSERT_VALUES); CHKERRQ(ierr);
            }
          }
        }
        localIdx++;
      }
    }
  }

  ierr = MatAssemblyBegin(ET, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(ET, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  ierr = PetscObjectViewFromOptions((PetscObject) ET, NULL, "-ET_mat_view"); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // generateET
