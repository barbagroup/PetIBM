#include "TairaColoniusSolver.h"
#include <petscsys.h>
#include <petscdmcomposite.h>
#include <iostream>
#include <string>

#include "TairaColonius/generateBNQ.inl"
#include "TairaColonius/generateR2.inl"
#include "TairaColonius/initialiseBodies.inl"

template <>
void TairaColoniusSolver<2>::initialise()
{
	PetscErrorCode ierr;
	PetscInt       rank, numProcs;
	PetscInt       m, n;
	const PetscInt *lxp, *lyp;

	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	MPI_Comm_size(PETSC_COMM_WORLD, &numProcs);

	initialiseBodies();
	createDMs();

	ierr = DMDAGetOwnershipRanges(pda, &lxp, &lyp, NULL); CHKERRV(ierr);
	ierr = DMDAGetInfo(pda, NULL, NULL, NULL, NULL, &m, &n, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRV(ierr);

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
	ierr = DMCompositeAddDM(lambdaPack, bda); CHKERRV(ierr);

	createVecs();

	PetscInt lambdaStart, lambdaEnd;

	ierr = VecGetOwnershipRange(lambda, &lambdaStart, &lambdaEnd); CHKERRV(ierr);
	ierr = DMDAGetCorners(pda, NULL, NULL, NULL, &m, &n, NULL); CHKERRV(ierr);
	startGlobalIndex = lambdaStart + m*n;

	MPI_Barrier(PETSC_COMM_WORLD);
	MPI_Allgather(&startGlobalIndex, 1, MPIU_INT, &startGlobalIndices.front(), 1, MPIU_INT, PETSC_COMM_WORLD);

	PetscInt globalIndex;

	for(PetscInt j=0; j<numProcs; j++)
	{
		globalIndex = startGlobalIndices[j];
		for(auto i=boundaryPointIndices[j].begin(); i!=boundaryPointIndices[j].end(); i++)
		{
			bodyGlobalIndices[*i] = globalIndex;
			globalIndex++;
		}
	}

	createLocalToGlobalMappingsFluxes();
	createLocalToGlobalMappingsLambda();
	initialiseMeshSpacings();
	initialiseFluxes();
	updateBoundaryGhosts();

	generateDiagonalMatrices();
	generateA();
	generateBNQ();
	generateQTBNQ();
	createKSPs();
}

template <>
void TairaColoniusSolver<3>::initialise()
{
	PetscErrorCode ierr;
	PetscInt       rank, numProcs;
	PetscInt       m, n, p;
	const PetscInt *lxp, *lyp, *lzp;

	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	MPI_Comm_size(PETSC_COMM_WORLD, &numProcs);

	initialiseBodies();
	createDMs();
	
	ierr = DMDAGetOwnershipRanges(pda, &lxp, &lyp, &lzp); CHKERRV(ierr);
	ierr = DMDAGetInfo(pda, NULL, NULL, NULL, NULL, &m, &n, &p, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRV(ierr);

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
	ierr = DMCompositeAddDM(lambdaPack, bda); CHKERRV(ierr);

	createVecs();

	PetscInt lambdaStart, lambdaEnd;

	ierr = VecGetOwnershipRange(lambda, &lambdaStart, &lambdaEnd); CHKERRV(ierr);
	ierr = DMDAGetCorners(pda, NULL, NULL, NULL, &m, &n, &p); CHKERRV(ierr);
	startGlobalIndex = lambdaStart + m*n*p;

	MPI_Barrier(PETSC_COMM_WORLD);
	MPI_Allgather(&startGlobalIndex, 1, MPIU_INT, &startGlobalIndices.front(), 1, MPIU_INT, PETSC_COMM_WORLD);

	PetscInt globalIndex;

	for(PetscInt j=0; j<numProcs; j++)
	{
		globalIndex = startGlobalIndices[j];
		for(auto i=boundaryPointIndices[j].begin(); i!=boundaryPointIndices[j].end(); i++)
		{
			bodyGlobalIndices[*i] = globalIndex;
			globalIndex++;
		}
	}

	createLocalToGlobalMappingsFluxes();
	createLocalToGlobalMappingsLambda();
	initialiseMeshSpacings();
	initialiseFluxes();
	updateBoundaryGhosts();

	generateDiagonalMatrices();
	generateA();
	generateBNQ();
	generateQTBNQ();
	createKSPs();
}

template <PetscInt dim>
void TairaColoniusSolver<dim>::finalise()
{
	NavierStokesSolver<dim>::finalise();
}

template class TairaColoniusSolver<2>;
template class TairaColoniusSolver<3>;