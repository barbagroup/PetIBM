#include "TairaColoniusSolver.h"
#include <petscsys.h>
#include <petscdmcomposite.h>
#include <iostream>
#include <string>

#include "TairaColonius/generateBNQ.inl"

template <>
void TairaColoniusSolver<2>::initialise()
{
	PetscErrorCode ierr;
	PetscInt       rank, numProcs;
	PetscInt       m, n;
	const PetscInt *lxp, *lyp;
	PetscInt       totalPoints = 51;
	
	h = 1.0/32;

	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	MPI_Comm_size(PETSC_COMM_WORLD, &numProcs);
	PetscPrintf(PETSC_COMM_WORLD, "rank: %d, numProcs: %d\n", rank, numProcs);

	numBoundaryPointsOnProcess.resize(numProcs);

	for(PetscInt l=0; l < totalPoints; l++)
	{
		x.push_back(0.5+0.25*cos(2*PETSC_PI*l/totalPoints));
		y.push_back(0.5+0.25*sin(2*PETSC_PI*l/totalPoints));
	}
	bodyGlobalIndices.resize(x.size());

	createDMs();

	ierr = DMDAGetOwnershipRanges(pda, &lxp, &lyp, NULL); CHKERRV(ierr);
	ierr = DMDAGetInfo(pda, NULL, NULL, NULL, NULL, &m, &n, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRV(ierr);

	boundaryPointIndices.resize(numProcs);
	startGlobalIndices.resize(numProcs);

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

	PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[%d]\n", rank);

	PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Body indices on each process\n", rank);
	for(PetscInt j=0; j<numProcs; j++)
	{
		PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[%d] ", j);
		for(auto i=boundaryPointIndices[j].begin(); i!=boundaryPointIndices[j].end(); i++)
		{
			PetscSynchronizedPrintf(PETSC_COMM_WORLD, "%d, ", *i);
		}
		PetscSynchronizedPrintf(PETSC_COMM_WORLD, "\n");
	}

	PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Number of body points on each process\n");
	for(PetscInt i=0; i<numProcs; i++)
		PetscSynchronizedPrintf(PETSC_COMM_WORLD, "%d, ", numBoundaryPointsOnProcess[i]);
	PetscSynchronizedPrintf(PETSC_COMM_WORLD, "\n");

	ierr = DMCreateGlobalVector(lambdaPack, &lambda); CHKERRV(ierr);

	PetscInt lambdaStart, lambdaEnd;

	ierr = VecGetOwnershipRange(lambda, &lambdaStart, &lambdaEnd); CHKERRV(ierr);
	ierr = DMDAGetCorners(pda, NULL, NULL, NULL, &m, &n, NULL); CHKERRV(ierr);
	startGlobalIndex = lambdaStart + m*n;

	PetscSynchronizedPrintf(PETSC_COMM_WORLD, "startGlobalIndex: %d\n", startGlobalIndex);	

	MPI_Barrier(PETSC_COMM_WORLD);
	MPI_Allgather(&startGlobalIndex, 1, MPIU_INT, &startGlobalIndices.front(), 1, MPIU_INT, PETSC_COMM_WORLD);

	PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Global indices of body points on each process\n");
	for(PetscInt i=0; i<numProcs; i++)
		PetscSynchronizedPrintf(PETSC_COMM_WORLD, "%d, ", startGlobalIndices[i]);
	PetscSynchronizedPrintf(PETSC_COMM_WORLD, "\n");

	PetscSynchronizedPrintf(PETSC_COMM_WORLD, "\n");
	PetscSynchronizedFlush(PETSC_COMM_WORLD);

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

	PetscPrintf(PETSC_COMM_WORLD, "Mapping\n");
	for(size_t i=0; i<bodyGlobalIndices.size(); i++)
	{
		PetscPrintf(PETSC_COMM_WORLD, "%d : %d\n", i, bodyGlobalIndices[i]);
	}

	createVecs();
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
}

template <PetscInt dim>
void TairaColoniusSolver<dim>::finalise()
{
	NavierStokesSolver<dim>::finalise();
}

template class TairaColoniusSolver<2>;
template class TairaColoniusSolver<3>;