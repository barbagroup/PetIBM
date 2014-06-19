template <>
PetscErrorCode TairaColoniusSolver<2>::generateET()
{
	PetscErrorCode ierr;
	PetscInt       numProcs;
	PetscInt       mstart, nstart, m, n;
	PetscInt       *d_nnz, *o_nnz;
	PetscInt       qStart, qEnd, qLocalSize, fStart, fEnd, fLocalSize;
	PetscInt       localIdx;
	PetscInt       row, col, value;
	PetscReal      xCoord, yCoord, h;
	PetscReal      disp[2];
	Vec            fGlobal;
	PetscLogEvent  GENERATE_ET;
	
	ierr = PetscLogEventRegister("generateET", 0, &GENERATE_ET); CHKERRQ(ierr);
	ierr = PetscLogEventBegin(GENERATE_ET, 0, 0, 0, 0); CHKERRQ(ierr);
	
	ierr = MPI_Comm_size(PETSC_COMM_WORLD, &numProcs); CHKERRQ(ierr);

	// ownership range of q
	ierr = VecGetOwnershipRange(q, &qStart, &qEnd); CHKERRQ(ierr);
	qLocalSize = qEnd-qStart;

	ierr = DMCompositeGetAccess(lambdaPack, lambda,  NULL, &fGlobal); CHKERRQ(ierr);

	// ownership range of f
	ierr = VecGetOwnershipRange(fGlobal, &fStart, &fEnd); CHKERRQ(ierr);
	fLocalSize = fEnd-fStart;

	// create arrays to store nnz values
	ierr = PetscMalloc(qLocalSize*sizeof(PetscInt), &d_nnz); CHKERRQ(ierr);
	ierr = PetscMalloc(qLocalSize*sizeof(PetscInt), &o_nnz); CHKERRQ(ierr);

	// determine the number of non-zeros in each row
	// in the diagonal and off-diagonal portions of the matrix
	localIdx = 0;
	// U
	ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
	for(PetscInt j=nstart; j<nstart+n; j++)
	{
		yCoord = 0.5*(mesh->y[j] + mesh->y[j+1]);
		for(PetscInt i=mstart; i<mstart+m; i++)
		{
			h = mesh->dx[i];
			xCoord = mesh->x[i+1];
			d_nnz[localIdx] = 0;
			o_nnz[localIdx] = 0;
			// ET portion
			PetscInt numPhi = 0;
			for(PetscInt procIdx=0; procIdx<numProcs; procIdx++)
			{
				numPhi += numPhiOnProcess[procIdx];
				for(auto l=boundaryPointIndices[procIdx].begin(); l!=boundaryPointIndices[procIdx].end(); l++)
				{
					if(isInfluenced(xCoord, yCoord, x[*l], y[*l], 1.5*h, disp))
					{
						col = globalIndexMapping[*l] - numPhi;
						(col>=fStart && col<fEnd)? d_nnz[localIdx]++ : o_nnz[localIdx]++;
					}
				}
			}
			localIdx++;
		}
	}
	// V
	ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
	for(PetscInt j=nstart; j<nstart+n; j++)
	{
		h = mesh->dy[j];
		yCoord = mesh->y[j+1];
		for(PetscInt i=mstart; i<mstart+m; i++)
		{
			xCoord = 0.5*(mesh->x[i] + mesh->x[i+1]);
			d_nnz[localIdx] = 0;
			o_nnz[localIdx] = 0;
			// ET portion
			PetscInt numPhi = 0;
			for(PetscInt procIdx=0; procIdx<numProcs; procIdx++)
			{
				numPhi += numPhiOnProcess[procIdx];
				for(auto l=boundaryPointIndices[procIdx].begin(); l!=boundaryPointIndices[procIdx].end(); l++)
				{
					if(isInfluenced(xCoord, yCoord, x[*l], y[*l], 1.5*h, disp))
					{
						col = globalIndexMapping[*l] - numPhi + 1;
						(col>=fStart && col<fEnd)? d_nnz[localIdx]++ : o_nnz[localIdx]++;
					}
				}
			}
			localIdx++;
		}
	}
	
	// allocate memory for the matrix
	ierr = MatCreate(PETSC_COMM_WORLD, &ET); CHKERRQ(ierr);
	ierr = MatSetType(ET, MATMPIAIJ); CHKERRQ(ierr);
	ierr = MatSetSizes(ET, qLocalSize, fLocalSize, PETSC_DETERMINE, PETSC_DETERMINE); CHKERRQ(ierr);
	ierr = MatMPIAIJSetPreallocation(ET, 0, d_nnz, 0, o_nnz); CHKERRQ(ierr);

	// deallocate d_nnz and o_nnz
	ierr = PetscFree(d_nnz); CHKERRQ(ierr);
	ierr = PetscFree(o_nnz); CHKERRQ(ierr);

	// assemble matrix Q
	localIdx = 0;
	// U
	ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
	for(PetscInt j=nstart; j<nstart+n; j++)
	{
		yCoord = 0.5*(mesh->y[j] + mesh->y[j+1]);
		for(PetscInt i=mstart; i<mstart+m; i++)
		{
			h = mesh->dx[i];
			xCoord = mesh->x[i+1];
			row = localIdx + qStart;
			// ET portion
			PetscInt numPhi = 0;
			for(PetscInt procIdx=0; procIdx<numProcs; procIdx++)
			{
				numPhi += numPhiOnProcess[procIdx];
				for(auto l=boundaryPointIndices[procIdx].begin(); l!=boundaryPointIndices[procIdx].end(); l++)
				{
					if(isInfluenced(xCoord, yCoord, x[*l], y[*l], 1.5*h, disp))
					{
						col  = globalIndexMapping[*l] - numPhi;
						value= h*delta(disp[0], disp[1], h);
						ierr = MatSetValue(ET, row, col, value, INSERT_VALUES); CHKERRQ(ierr);
					}
				}
			}
			localIdx++;
		}
	}
	// V
	ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
	for(PetscInt j=nstart; j<nstart+n; j++)
	{
		h = mesh->dy[j];
		yCoord = mesh->y[j+1];
		for(PetscInt i=mstart; i<mstart+m; i++)
		{
			xCoord = 0.5*(mesh->x[i] + mesh->x[i+1]);
			row = localIdx + qStart;
			// ET portion
			PetscInt numPhi = 0;
			for(PetscInt procIdx=0; procIdx<numProcs; procIdx++)
			{
				numPhi += numPhiOnProcess[procIdx];
				for(auto l=boundaryPointIndices[procIdx].begin(); l!=boundaryPointIndices[procIdx].end(); l++)
				{
					if(isInfluenced(xCoord, yCoord, x[*l], y[*l], 1.5*h, disp))
					{
						col = globalIndexMapping[*l] - numPhi + 1;
						value= h*delta(disp[0], disp[1], h);
						ierr = MatSetValue(ET, row, col, value, INSERT_VALUES); CHKERRQ(ierr);
					}
				}
			}
			localIdx++;
		}
	}

	ierr = MatAssemblyBegin(ET, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(ET, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	
	ierr = PetscLogEventEnd(GENERATE_ET, 0, 0, 0, 0); CHKERRQ(ierr);

	return 0;
}

template <>
PetscErrorCode TairaColoniusSolver<3>::generateET()
{
	PetscErrorCode ierr;
	PetscInt       numProcs;
	PetscInt       mstart, nstart, pstart, m, n, p;
	PetscInt       *d_nnz, *o_nnz;
	PetscInt       qStart, qEnd, qLocalSize, fStart, fEnd, fLocalSize;
	PetscInt       localIdx;
	PetscInt       row, col, value;
	PetscReal      xCoord, yCoord, zCoord, h;
	PetscReal      disp[3];
	Vec            fGlobal;
	PetscLogEvent  GENERATE_ET;
	
	ierr = PetscLogEventRegister("generateET", 0, &GENERATE_ET); CHKERRQ(ierr);
	ierr = PetscLogEventBegin(GENERATE_ET, 0, 0, 0, 0); CHKERRQ(ierr);

	ierr = MPI_Comm_size(PETSC_COMM_WORLD, &numProcs); CHKERRQ(ierr);

	// ownership range of q
	ierr = VecGetOwnershipRange(q, &qStart, &qEnd); CHKERRQ(ierr);
	qLocalSize = qEnd-qStart;

	ierr = DMCompositeGetAccess(lambdaPack, lambda, NULL, &fGlobal); CHKERRQ(ierr);

	// ownership range of f
	ierr = VecGetOwnershipRange(fGlobal, &fStart, &fEnd); CHKERRQ(ierr);
	fLocalSize = fEnd-fStart;

	// create arrays to store nnz values
	ierr = PetscMalloc(qLocalSize*sizeof(PetscInt), &d_nnz); CHKERRQ(ierr);
	ierr = PetscMalloc(qLocalSize*sizeof(PetscInt), &o_nnz); CHKERRQ(ierr);

	// determine the number of non-zeros in each row
	// in the diagonal and off-diagonal portions of the matrix
	localIdx = 0;
	// U
	ierr = DMDAGetCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
	for(PetscInt k=pstart; k<pstart+p; k++)
	{
		zCoord = 0.5*(mesh->z[k] + mesh->z[k+1]);
		for(PetscInt j=nstart; j<nstart+n; j++)
		{
			yCoord = 0.5*(mesh->y[j] + mesh->y[j+1]);
			for(PetscInt i=mstart; i<mstart+m; i++)
			{
				h = mesh->dx[i];
				xCoord = mesh->x[i+1];
				d_nnz[localIdx] = 0;
				o_nnz[localIdx] = 0;
				// ET portion
				PetscInt numPhi = 0;
				for(PetscInt procIdx=0; procIdx<numProcs; procIdx++)
				{
					numPhi += numPhiOnProcess[procIdx];
					for(auto l=boundaryPointIndices[procIdx].begin(); l!=boundaryPointIndices[procIdx].end(); l++)
					{
						if(isInfluenced(xCoord, yCoord, zCoord, x[*l], y[*l], z[*l], 1.5*h, disp))
						{
							col = globalIndexMapping[*l] - numPhi;
							(col>=fStart && col<fEnd)? d_nnz[localIdx]++ : o_nnz[localIdx]++;
						}
					}
				}
				localIdx++;
			}
		}
	}
	// V
	ierr = DMDAGetCorners(vda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
	for(PetscInt k=pstart; k<pstart+p; k++)
	{
		zCoord = 0.5*(mesh->z[k] + mesh->z[k+1]);
		for(PetscInt j=nstart; j<nstart+n; j++)
		{
			h = mesh->dy[j];
			yCoord = mesh->y[j+1];
			for(PetscInt i=mstart; i<mstart+m; i++)
			{
				xCoord = 0.5*(mesh->x[i] + mesh->x[i+1]);
				d_nnz[localIdx] = 0;
				o_nnz[localIdx] = 0;
				// ET portion
				PetscInt numPhi = 0;
				for(PetscInt procIdx=0; procIdx<numProcs; procIdx++)
				{
					numPhi += numPhiOnProcess[procIdx];
					for(auto l=boundaryPointIndices[procIdx].begin(); l!=boundaryPointIndices[procIdx].end(); l++)
					{
						if(isInfluenced(xCoord, yCoord, zCoord, x[*l], y[*l], z[*l], 1.5*h, disp))
						{
							col = globalIndexMapping[*l] - numPhi + 1;
							(col>=fStart && col<fEnd)? d_nnz[localIdx]++ : o_nnz[localIdx]++;
						}
					}
				}
				localIdx++;
			}
		}
	}
	// W
	ierr = DMDAGetCorners(wda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
	for(PetscInt k=pstart; k<pstart+p; k++)
	{
		h = mesh->dz[k];
		zCoord = mesh->z[k+1];
		for(PetscInt j=nstart; j<nstart+n; j++)
		{
			yCoord = 0.5*(mesh->y[j] + mesh->y[j+1]);
			for(PetscInt i=mstart; i<mstart+m; i++)
			{
				xCoord = 0.5*(mesh->x[i] + mesh->x[i+1]);
				d_nnz[localIdx] = 0;
				o_nnz[localIdx] = 0;
				// ET portion
				PetscInt numPhi = 0;
				for(PetscInt procIdx=0; procIdx<numProcs; procIdx++)
				{
					numPhi += numPhiOnProcess[procIdx];
					for(auto l=boundaryPointIndices[procIdx].begin(); l!=boundaryPointIndices[procIdx].end(); l++)
					{
						if(isInfluenced(xCoord, yCoord, zCoord, x[*l], y[*l], z[*l], 1.5*h, disp))
						{
							col = globalIndexMapping[*l] - numPhi + 2;
							(col>=fStart && col<fEnd)? d_nnz[localIdx]++ : o_nnz[localIdx]++;
						}
					}
				}
				localIdx++;
			}
		}
	}
	
	// allocate memory for the matrix
	ierr = MatCreate(PETSC_COMM_WORLD, &ET); CHKERRQ(ierr);
	ierr = MatSetType(ET, MATMPIAIJ); CHKERRQ(ierr);
	ierr = MatSetSizes(ET, qLocalSize, fLocalSize, PETSC_DETERMINE, PETSC_DETERMINE); CHKERRQ(ierr);
	ierr = MatMPIAIJSetPreallocation(ET, 0, d_nnz, 0, o_nnz); CHKERRQ(ierr);

	// deallocate d_nnz and o_nnz
	ierr = PetscFree(d_nnz); CHKERRQ(ierr);
	ierr = PetscFree(o_nnz); CHKERRQ(ierr);

	// assemble matrix Q
	localIdx = 0;
	// U
	ierr = DMDAGetCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
	for(PetscInt k=pstart; k<pstart+p; k++)
	{
		zCoord = 0.5*(mesh->z[k] + mesh->z[k+1]);
		for(PetscInt j=nstart; j<nstart+n; j++)
		{
			yCoord = 0.5*(mesh->y[j] + mesh->y[j+1]);
			for(PetscInt i=mstart; i<mstart+m; i++)
			{
				h = mesh->dx[i];
				xCoord = mesh->x[i+1];
				row = localIdx + qStart;
				// ET portion
				PetscInt numPhi = 0;
				for(PetscInt procIdx=0; procIdx<numProcs; procIdx++)
				{
					numPhi += numPhiOnProcess[procIdx];
					for(auto l=boundaryPointIndices[procIdx].begin(); l!=boundaryPointIndices[procIdx].end(); l++)
					{
						if(isInfluenced(xCoord, yCoord, zCoord, x[*l], y[*l], z[*l], 1.5*h, disp))
						{
							col  = globalIndexMapping[*l] - numPhi;
							value= h*delta(disp[0], disp[1], disp[2], h);
							ierr = MatSetValue(ET, row, col, value, INSERT_VALUES); CHKERRQ(ierr);
						}
					}
				}
				localIdx++;
			}
		}
	}
	// V
	ierr = DMDAGetCorners(vda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
	for(PetscInt k=pstart; k<pstart+p; k++)
	{
		zCoord = 0.5*(mesh->z[k] + mesh->z[k+1]);
		for(PetscInt j=nstart; j<nstart+n; j++)
		{
			h = mesh->dy[j];
			yCoord = mesh->y[j+1];
			for(PetscInt i=mstart; i<mstart+m; i++)
			{
				xCoord = 0.5*(mesh->x[i] + mesh->x[i+1]);
				row = localIdx + qStart;
				// ET portion
				PetscInt numPhi = 0;
				for(PetscInt procIdx=0; procIdx<numProcs; procIdx++)
				{
					numPhi += numPhiOnProcess[procIdx];
					for(auto l=boundaryPointIndices[procIdx].begin(); l!=boundaryPointIndices[procIdx].end(); l++)
					{
						if(isInfluenced(xCoord, yCoord, zCoord, x[*l], y[*l], z[*l], 1.5*h, disp))
						{
							col  = globalIndexMapping[*l] - numPhi + 1;
							value= h*delta(disp[0], disp[1], disp[2], h);
							ierr = MatSetValue(ET, row, col, value, INSERT_VALUES); CHKERRQ(ierr);
						}
					}
				}
				localIdx++;
			}
		}
	}
	// W
	ierr = DMDAGetCorners(wda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
	for(PetscInt k=pstart; k<pstart+p; k++)
	{
		h = mesh->dz[k];
		zCoord = mesh->z[k+1];
		for(PetscInt j=nstart; j<nstart+n; j++)
		{
			yCoord = 0.5*(mesh->y[j] + mesh->y[j+1]);
			for(PetscInt i=mstart; i<mstart+m; i++)
			{
				xCoord = 0.5*(mesh->x[i] + mesh->x[i+1]);
				row = localIdx + qStart;
				// ET portion
				PetscInt numPhi = 0;
				for(PetscInt procIdx=0; procIdx<numProcs; procIdx++)
				{
					numPhi += numPhiOnProcess[procIdx];
					for(auto l=boundaryPointIndices[procIdx].begin(); l!=boundaryPointIndices[procIdx].end(); l++)
					{
						if(isInfluenced(xCoord, yCoord, zCoord, x[*l], y[*l], z[*l], 1.5*h, disp))
						{
							col  = globalIndexMapping[*l] - numPhi + 2;
							value= h*delta(disp[0], disp[1], disp[2], h);
							ierr = MatSetValue(ET, row, col, value, INSERT_VALUES); CHKERRQ(ierr);
						}
					}
				}
				localIdx++;
			}
		}
	}

	ierr = MatAssemblyBegin(ET, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(ET, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	
	ierr = PetscLogEventEnd(GENERATE_ET, 0, 0, 0, 0); CHKERRQ(ierr);

	return 0;
}