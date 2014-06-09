template <>
PetscErrorCode TairaColoniusSolver<2>::createDMs()
{
	NavierStokesSolver<2>::createDMs();

	PetscErrorCode ierr;
	PetscInt       m, n;
	const PetscInt *lxp, *lyp;

	ierr = DMDAGetOwnershipRanges(pda, &lxp, &lyp, NULL); CHKERRQ(ierr);
	ierr = DMDAGetInfo(pda, NULL, NULL, NULL, NULL, &m, &n, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);

	PetscInt xStart, yStart, xEnd, yEnd,
	         procIdx = 0;

	yStart = 0;
	for(PetscInt j=0; j<n; j++)
	{
		yEnd   = yStart + lyp[j];
		xStart = 0;
		for(PetscInt i=0; i<m; i++)
		{
			procIdx = j*m + i;
			xEnd = xStart + lxp[i];
			for(size_t l=0; l<x.size(); l++)
			{
				if(x[l]>=mesh->x[xStart] && x[l]<mesh->x[xEnd] && y[l]>=mesh->y[yStart] && y[l]<mesh->y[yEnd])
				{
					boundaryPointIndices[procIdx].push_back(l);
					numBoundaryPointsOnProcess[procIdx]++;
				}
			}
			xStart = xEnd;
		}
		yStart = yEnd;
	}

	ierr = DMDACreate1d(PETSC_COMM_WORLD, DMDA_BOUNDARY_NONE, x.size(), 2, 0, &numBoundaryPointsOnProcess.front(), &bda);
	ierr = DMCompositeAddDM(lambdaPack, bda); CHKERRQ(ierr);

	return 0;
}

template <>
PetscErrorCode TairaColoniusSolver<3>::createDMs()
{
	NavierStokesSolver<3>::createDMs();

	PetscErrorCode ierr;
	PetscInt       m, n, p;
	const PetscInt *lxp, *lyp, *lzp;

	ierr = DMDAGetOwnershipRanges(pda, &lxp, &lyp, &lzp); CHKERRQ(ierr);
	ierr = DMDAGetInfo(pda, NULL, NULL, NULL, NULL, &m, &n, &p, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);

	PetscInt xStart, yStart, zStart, xEnd, yEnd, zEnd,
	         procIdx = 0;
	
	zStart = 0;
	for(PetscInt k=0; k<p; k++)
	{
		zEnd = zStart + lzp[k];
		yStart = 0;
		for(PetscInt j=0; j<n; j++)
		{
			yEnd   = yStart + lyp[j];
			xStart = 0;
			for(PetscInt i=0; i<m; i++)
			{
				procIdx = k*m*n + j*m + i;
				xEnd = xStart + lxp[i];
				for(size_t l=0; l<x.size(); l++)
				{
					if(x[l]>=mesh->x[xStart] && x[l]<mesh->x[xEnd] && y[l]>=mesh->y[yStart] && y[l]<mesh->y[yEnd] && z[l]>=mesh->z[zStart] && z[l]<mesh->z[zEnd])
					{
						boundaryPointIndices[procIdx].push_back(l);
						numBoundaryPointsOnProcess[procIdx]++;
					}
				}
				xStart = xEnd;
			}
			yStart = yEnd;
		}
		zStart = zEnd;
	}
	
	ierr = DMDACreate1d(PETSC_COMM_WORLD, DMDA_BOUNDARY_NONE, x.size(), 3, 0, &numBoundaryPointsOnProcess.front(), &bda);
	ierr = DMCompositeAddDM(lambdaPack, bda); CHKERRQ(ierr);

	return 0;
}