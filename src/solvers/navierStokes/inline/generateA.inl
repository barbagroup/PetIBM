/***************************************************************************//**
 * \file generateA.inl
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the method `generateA` of the class `NavierStokesSolver`.
 */


void getColumnIndices(PetscReal **mappingLocalToGlobal, PetscInt i, PetscInt j, PetscInt *cols)
{
  cols[0] = mappingLocalToGlobal[j][i];
  cols[1] = mappingLocalToGlobal[j][i-1];
  cols[2] = mappingLocalToGlobal[j][i+1];
  cols[3] = mappingLocalToGlobal[j-1][i];  
  cols[4] = mappingLocalToGlobal[j+1][i];
} // getColumns

void getColumnIndices(PetscReal ***mappingLocalToGlobal, PetscInt i, PetscInt j, PetscInt k, PetscInt *cols)
{
  cols[0] = mappingLocalToGlobal[k][j][i];
  cols[1] = mappingLocalToGlobal[k][j][i-1];
  cols[2] = mappingLocalToGlobal[k][j][i+1];
  cols[3] = mappingLocalToGlobal[k][j-1][i]; 
  cols[4] = mappingLocalToGlobal[k][j+1][i];
  cols[5] = mappingLocalToGlobal[k-1][j][i];
  cols[6] = mappingLocalToGlobal[k+1][j][i];
} // getColumns

void getCoefficients(PetscReal dxMinus, PetscReal dxPlus, PetscReal dyMinus, PetscReal dyPlus, PetscReal *values)
{
  values[0] = -(2.0/dxMinus/dxPlus + 2.0/dyMinus/dyPlus);
  values[1] = 2.0/dxMinus/(dxMinus + dxPlus);
  values[2] = 2.0/ dxPlus/(dxMinus + dxPlus);
  values[3] = 2.0/dyMinus/(dyMinus + dyPlus);
  values[4] = 2.0/ dyPlus/(dyMinus + dyPlus);
} // getCoefficients

void getCoefficients(PetscReal dxMinus, PetscReal dxPlus, PetscReal dyMinus, PetscReal dyPlus, PetscReal dzMinus, PetscReal dzPlus, PetscReal *values)
{
  getCoefficients(dxMinus, dxPlus, dyMinus, dyPlus, values);  
  values[0] += (-2.0/dzMinus/dzPlus);
  values[5] = 2.0/dzMinus/(dzMinus + dzPlus);
  values[6] = 2.0/ dzPlus/(dzMinus + dzPlus);
} // getCoefficients


/**
 * \brief Assembles the matrix that results from implicit contributions of the 
 *        dicretized momentum equations.
 *
 * The matix is composed of the implicit coefficients from the time-derivative
 * as well as the implicit coefficients from the diffusive terms.
 * Moreover, the matrix a diagonally scaled by the matrices \f$ \hat{M} \f$ 
 * and \f$ R^{-1} \f$.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::generateA()
{
  return 0;
} // generateA


// two-dimensional specialization
template <>
PetscErrorCode NavierStokesSolver<2>::generateA()
{
  PetscErrorCode ierr;

  PetscInt i, j,           // loop indices
           M, N,           // global number of nodes along each direction
           m, n,           // local number of nodes along each direction
           mstart, nstart; // starting indices
  
  // get local size of the fluxes
  PetscInt qStart, qEnd, qLocalSize;
  ierr = VecGetOwnershipRange(q, &qStart, &qEnd); CHKERRQ(ierr);
  qLocalSize = qEnd-qStart;

  // create arrays to store number of non-zeros (nnz) values
  PetscInt *d_nnz, // nnz on diagonal
           *o_nnz; // nnz off diagonal
  ierr = PetscMalloc(qLocalSize*sizeof(PetscInt), &d_nnz); CHKERRQ(ierr);
  ierr = PetscMalloc(qLocalSize*sizeof(PetscInt), &o_nnz); CHKERRQ(ierr);

  // determine the number of non-zeros row by row
  PetscInt localIdx = 0;
  PetscInt cols[5];
  PetscReal values[5];
  // rows corresponding to fluxes in x-direction
  PetscReal **uMappingArray;
  ierr = DMDAVecGetArray(uda, uMapping, &uMappingArray); CHKERRQ(ierr);
  ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
  for (j=nstart; j<nstart+n; j++)
  {
    for (i=mstart; i<mstart+m; i++)
    {
      getColumnIndices(uMappingArray, i, j, cols);
      countNumNonZeros(cols, 5, qStart, qEnd, d_nnz[localIdx], o_nnz[localIdx]);
      localIdx++;
    }
  }
  ierr = DMDAVecRestoreArray(uda, uMapping, &uMappingArray); CHKERRQ(ierr);
  // rows corresponding to fluxes in y-direction
  PetscReal **vMappingArray;
  ierr = DMDAVecGetArray(vda, vMapping, &vMappingArray); CHKERRQ(ierr);
  ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
  for (j=nstart; j<nstart+n; j++)
  {
    for (i=mstart; i<mstart+m; i++)
    {
      getColumnIndices(vMappingArray, i, j, cols);
      countNumNonZeros(cols, 5, qStart, qEnd, d_nnz[localIdx], o_nnz[localIdx]);
      localIdx++;
    }
  }
  ierr = DMDAVecRestoreArray(vda, vMapping, &vMappingArray); CHKERRQ(ierr);

  // create and allocate memory for matrix A
  ierr = MatCreate(PETSC_COMM_WORLD, &A); CHKERRQ(ierr);
  ierr = MatSetSizes(A, qLocalSize, qLocalSize, PETSC_DETERMINE, PETSC_DETERMINE); CHKERRQ(ierr);
  ierr = MatSetFromOptions(A); CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(A, 0, d_nnz); CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(A, 0, d_nnz, 0, o_nnz); CHKERRQ(ierr);

  // deallocate d_nnz and o_nnz
  ierr = PetscFree(d_nnz); CHKERRQ(ierr);
  ierr = PetscFree(o_nnz); CHKERRQ(ierr);

  // assemble matrix A row by row
  PetscReal dxMinus, dxPlus, dyMinus, dyPlus;
  // rows corresponding to fluxes in x-direction
  ierr = DMDAVecGetArray(uda, uMapping, &uMappingArray); CHKERRQ(ierr);
  ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
  ierr = DMDAGetInfo(uda, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
  for (j=nstart; j<nstart+n; j++)
  {
    dyMinus = (j == 0) ? 0.5*mesh->dy[0] : 0.5*(mesh->dy[j-1]+mesh->dy[j]);
    dyPlus = (j == N-1) ? 0.5*mesh->dy[N-1] : 0.5*(mesh->dy[j]+mesh->dy[j+1]);
    for (i=mstart; i<mstart+m; i++)
    {
      dxMinus = mesh->dx[i];
      dxPlus = mesh->dx[i+1];
      getColumnIndices(uMappingArray, i, j, cols);
      getCoefficients(dxMinus, dxPlus, dyMinus, dyPlus, values);
      ierr = MatSetValues(A, 1, &cols[0], 5, cols, values, INSERT_VALUES); CHKERRQ(ierr);
    }
  }
  ierr = DMDAVecRestoreArray(uda, uMapping, &uMappingArray); CHKERRQ(ierr);
  // row corresponding to fluxes in y-direction
  ierr = DMDAVecGetArray(vda, vMapping, &vMappingArray); CHKERRQ(ierr);
  ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
  ierr = DMDAGetInfo(vda, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
  for (j=nstart; j<nstart+n; j++)
  {
    dyMinus = mesh->dy[j];
    dyPlus = mesh->dy[j+1];
    for (i=mstart; i<mstart+m; i++)
    {
      dxMinus = (i == 0) ? 0.5*mesh->dx[0] : 0.5*(mesh->dx[i-1]+mesh->dx[i]);
      dxPlus = (i == M-1) ? 0.5*mesh->dx[M-1] : 0.5*(mesh->dx[i]+mesh->dx[i+1]);
      getColumnIndices(vMappingArray, i, j, cols);
      getCoefficients(dxMinus, dxPlus, dyMinus, dyPlus, values);
      ierr = MatSetValues(A, 1, &cols[0], 5, cols, values, INSERT_VALUES); CHKERRQ(ierr);
    }
  }
  ierr = DMDAVecRestoreArray(vda, vMapping, &vMappingArray); CHKERRQ(ierr);

  ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  PetscReal alpha = parameters->diffusion.coefficients[0]; // implicit diffusion coefficient
  ierr = MatScale(A, -alpha*flow->nu); CHKERRQ(ierr);
  ierr = MatShift(A, 1.0/parameters->dt); CHKERRQ(ierr);
  ierr = MatDiagonalScale(A, MHat, RInv);

  return 0;
} // generateA


// three-dimensional specialization
template <>
PetscErrorCode NavierStokesSolver<3>::generateA()
{
  PetscErrorCode ierr;

  PetscInt i, j, k,                // loop indices
           M, N, P,                // global number of nodes along each direction
           m, n, p,                // local number of nodes along each direction
           mstart, nstart, pstart; // startting indices

  // get local size of the fluxes
  PetscInt qStart, qEnd, qLocalSize;
  ierr = VecGetOwnershipRange(q, &qStart, &qEnd); CHKERRQ(ierr);
  qLocalSize = qEnd-qStart;

  // create arrays to store number of non-zero (nnz) values
  PetscInt *d_nnz, // nnz on diagonal
           *o_nnz; // nnz off diagonal
  ierr = PetscMalloc(qLocalSize*sizeof(PetscInt), &d_nnz); CHKERRQ(ierr);
  ierr = PetscMalloc(qLocalSize*sizeof(PetscInt), &o_nnz); CHKERRQ(ierr);

  // determine nnz row by row
  PetscInt localIdx = 0;
  PetscInt cols[7];
  PetscReal values[7];
  // rows corresponding to fluxes in x-direction
  PetscReal ***uMappingArray;
  ierr = DMDAVecGetArray(uda, uMapping, &uMappingArray); CHKERRQ(ierr);
  ierr = DMDAGetCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  for (k=pstart; k<pstart+p; k++)
  {
    for (j=nstart; j<nstart+n; j++)
    {
      for (i=mstart; i<mstart+m; i++)
      {
        getColumnIndices(uMappingArray, i, j, k, cols);
        countNumNonZeros(cols, 7, qStart, qEnd, d_nnz[localIdx], o_nnz[localIdx]);
        localIdx++;
      }
    }
  }
  ierr = DMDAVecRestoreArray(uda, uMapping, &uMappingArray); CHKERRQ(ierr);
  // rows corresponding to fluxes in y-direction
  PetscReal ***vMappingArray;
  ierr = DMDAVecGetArray(vda, vMapping, &vMappingArray); CHKERRQ(ierr);
  ierr = DMDAGetCorners(vda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  for (k=pstart; k<pstart+p; k++)
  {
    for (j=nstart; j<nstart+n; j++)
    {
      for (i=mstart; i<mstart+m; i++)
      {
        getColumnIndices(vMappingArray, i, j, k, cols);
        countNumNonZeros(cols, 7, qStart, qEnd, d_nnz[localIdx], o_nnz[localIdx]);
        localIdx++;
      }
    }
  }
  ierr = DMDAVecRestoreArray(vda, vMapping, &vMappingArray); CHKERRQ(ierr);
  // rows corresponding to fluxes in z-direction
  PetscReal ***wMappingArray;
  ierr = DMDAVecGetArray(wda, wMapping, &wMappingArray); CHKERRQ(ierr);
  ierr = DMDAGetCorners(wda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  for (k=pstart; k<pstart+p; k++)
  {
    for (j=nstart; j<nstart+n; j++)
    {
      for (i=mstart; i<mstart+m; i++)
      {
        getColumnIndices(wMappingArray, i, j, k, cols);
        countNumNonZeros(cols, 7, qStart, qEnd, d_nnz[localIdx], o_nnz[localIdx]);
        localIdx++;
      }
    }
  }
  ierr = DMDAVecRestoreArray(wda, wMapping, &wMappingArray); CHKERRQ(ierr);

  // create and allocate memory for matrix A
  ierr = MatCreate(PETSC_COMM_WORLD, &A); CHKERRQ(ierr);
  ierr = MatSetSizes(A, qLocalSize, qLocalSize, PETSC_DETERMINE, PETSC_DETERMINE); CHKERRQ(ierr);
  ierr = MatSetFromOptions(A); CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(A, 0, d_nnz); CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(A, 0, d_nnz, 0, o_nnz); CHKERRQ(ierr);

  // deallocate d_nnz and o_nnz
  ierr = PetscFree(d_nnz); CHKERRQ(ierr);
  ierr = PetscFree(o_nnz); CHKERRQ(ierr);

  // assemble matrix A row by row
  PetscReal dxMinus, dxPlus, dyMinus, dyPlus, dzMinus, dzPlus;
  // rows corresponding to fluxes in x-direction
  ierr = DMDAVecGetArray(uda, uMapping, &uMappingArray); CHKERRQ(ierr);
  ierr = DMDAGetCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  ierr = DMDAGetInfo(uda, NULL, &M, &N, &P, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
  for (k=pstart; k<pstart+p; k++)
  {
    dzMinus = (k == 0) ? 0.5*mesh->dz[0] : 0.5*(mesh->dz[k-1]+mesh->dz[k]);
    dzPlus = (k == P-1) ? 0.5*mesh->dz[P-1] : 0.5*(mesh->dz[k]+mesh->dz[k+1]);
    for (j=nstart; j<nstart+n; j++)
    {
      dyMinus = (j == 0) ? 0.5*mesh->dy[0] : 0.5*(mesh->dy[j-1]+mesh->dy[j]);
      dyPlus = (j == N-1) ? 0.5*mesh->dy[N-1] : 0.5*(mesh->dy[j]+mesh->dy[j+1]);
      for (i=mstart; i<mstart+m; i++)
      {
        dxMinus = mesh->dx[i];
        dxPlus = mesh->dx[i+1];
        getColumnIndices(uMappingArray, i, j, k, cols);
        getCoefficients(dxMinus, dxPlus, dyMinus, dyPlus, dzMinus, dzPlus, values);
        ierr = MatSetValues(A, 1, &cols[0], 7, cols, values, INSERT_VALUES); CHKERRQ(ierr);
      }
    }
  }
  ierr = DMDAVecRestoreArray(uda, uMapping, &uMappingArray); CHKERRQ(ierr);
  // rows corresponding to fluxes in y-direction
  ierr = DMDAVecGetArray(vda, vMapping, &vMappingArray); CHKERRQ(ierr);
  ierr = DMDAGetCorners(vda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  ierr = DMDAGetInfo(vda, NULL, &M, &N, &P, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
  for (k=pstart; k<pstart+p; k++)
  {
    dzMinus = (k == 0) ? 0.5*mesh->dz[0] : 0.5*(mesh->dz[k-1]+mesh->dz[k]);
    dzPlus = (k == P-1) ? 0.5*mesh->dz[P-1] : 0.5*(mesh->dz[k]+mesh->dz[k+1]);
    for (j=nstart; j<nstart+n; j++)
    {
      dyMinus = mesh->dy[j];
      dyPlus = mesh->dy[j+1];
      for (i=mstart; i<mstart+m; i++)
      {
        dxMinus = (i == 0) ? 0.5*mesh->dx[0] : 0.5*(mesh->dx[i-1]+mesh->dx[i]);
        dxPlus = (i == M-1) ? 0.5*mesh->dx[M-1] : 0.5*(mesh->dx[i]+mesh->dx[i+1]);
        getColumnIndices(vMappingArray, i, j, k, cols);
        getCoefficients(dxMinus, dxPlus, dyMinus, dyPlus, dzMinus, dzPlus, values);
        ierr = MatSetValues(A, 1, &cols[0], 7, cols, values, INSERT_VALUES); CHKERRQ(ierr);
      }
    }
  }
  ierr = DMDAVecRestoreArray(vda, vMapping, &vMappingArray); CHKERRQ(ierr);
  // rows corresponding to fluxes in z-direction
  ierr = DMDAVecGetArray(wda, wMapping, &wMappingArray); CHKERRQ(ierr);
  ierr = DMDAGetCorners(wda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
  ierr = DMDAGetInfo(wda, NULL, &M, &N, &P, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
  for (k=pstart; k<pstart+p; k++)
  {
    dzMinus = mesh->dz[k];
    dzPlus = mesh->dz[k+1];
    for (j=nstart; j<nstart+n; j++)
    {
      dyMinus = (j == 0) ? 0.5*mesh->dy[0] : 0.5*(mesh->dy[j-1]+mesh->dy[j]);
      dyPlus = (j == N-1) ? 0.5*mesh->dy[N-1] : 0.5*(mesh->dy[j]+mesh->dy[j+1]);
      for (i=mstart; i<mstart+m; i++)
      {
        dxMinus = (i == 0) ? 0.5*mesh->dx[0] : 0.5*(mesh->dx[i-1]+mesh->dx[i]);
        dxPlus = (i == M-1) ? 0.5*mesh->dx[M-1] : 0.5*(mesh->dx[i]+mesh->dx[i+1]);
        getColumnIndices(wMappingArray, i, j, k, cols);
        getCoefficients(dxMinus, dxPlus, dyMinus, dyPlus, dzMinus, dzPlus, values);
        ierr = MatSetValues(A, 1, &cols[0], 7, cols, values, INSERT_VALUES); CHKERRQ(ierr);
      }
    }
  }
  ierr = DMDAVecRestoreArray(wda, wMapping, &wMappingArray); CHKERRQ(ierr);

  ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  PetscReal alpha = parameters->diffusion.coefficients[0]; // implicit diffusion coefficient
  ierr = MatScale(A, -alpha*flow->nu); CHKERRQ(ierr);
  ierr = MatShift(A, 1.0/parameters->dt); CHKERRQ(ierr);
  ierr = MatDiagonalScale(A, MHat, RInv); CHKERRQ(ierr);

  return 0;
} // generateA