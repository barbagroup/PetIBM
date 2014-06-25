#include "NavierStokesSolver.h"
#include <petscdmcomposite.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <sys/stat.h>

template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::initialize()
{
	PetscErrorCode ierr;

	ierr = PetscLogStagePush(stageInitialize); CHKERRQ(ierr);
	ierr = createDMs(); CHKERRQ(ierr);
	ierr = initializeCommon(); CHKERRQ(ierr);
	ierr = PetscLogStagePop(); CHKERRQ(ierr);

	return 0;
}

template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::initializeCommon()
{
	PetscErrorCode ierr;

	ierr = createVecs(); CHKERRQ(ierr);
	
	initializeMeshSpacings();
	ierr = initializeFluxes(); CHKERRQ(ierr);
	ierr = initializeLambda(); CHKERRQ(ierr);
	ierr = updateBoundaryGhosts(); CHKERRQ(ierr);

	ierr = createLocalToGlobalMappingsFluxes(); CHKERRQ(ierr);
	ierr = createLocalToGlobalMappingsLambda(); CHKERRQ(ierr);

	ierr = generateDiagonalMatrices(); CHKERRQ(ierr);
	ierr = generateA(); CHKERRQ(ierr);
	ierr = generateBNQ(); CHKERRQ(ierr);
	ierr = generateQTBNQ(); CHKERRQ(ierr);
	ierr = createKSPs(); CHKERRQ(ierr);

	return 0;
}

template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::finalize()
{
	PetscErrorCode ierr;
	
	// DMs
	if(pda!=PETSC_NULL) {ierr = DMDestroy(&pda); CHKERRQ(ierr);}
	if(uda!=PETSC_NULL) {ierr = DMDestroy(&uda); CHKERRQ(ierr);}
	if(vda!=PETSC_NULL) {ierr = DMDestroy(&vda); CHKERRQ(ierr);}
	if(wda!=PETSC_NULL) {ierr = DMDestroy(&wda); CHKERRQ(ierr);}
	if(qPack!=PETSC_NULL){ierr = DMDestroy(&qPack); CHKERRQ(ierr);}
	if(lambdaPack!=PETSC_NULL){ierr = DMDestroy(&lambdaPack); CHKERRQ(ierr);}
	
	// Vecs
	if(q!=PETSC_NULL)    {ierr = VecDestroy(&q); CHKERRQ(ierr);}
	if(qStar!=PETSC_NULL){ierr = VecDestroy(&qStar); CHKERRQ(ierr);}
	
	if(qxLocal!=PETSC_NULL){ierr = VecDestroy(&qxLocal); CHKERRQ(ierr);}
	if(qyLocal!=PETSC_NULL){ierr = VecDestroy(&qyLocal); CHKERRQ(ierr);}
	if(qzLocal!=PETSC_NULL){ierr = VecDestroy(&qzLocal); CHKERRQ(ierr);}

	if(H!=PETSC_NULL)   {ierr = VecDestroy(&H); CHKERRQ(ierr);}
	if(rn!=PETSC_NULL)  {ierr = VecDestroy(&rn); CHKERRQ(ierr);}
	if(bc1!=PETSC_NULL) {ierr = VecDestroy(&bc1); CHKERRQ(ierr);}
	if(rhs1!=PETSC_NULL){ierr = VecDestroy(&rhs1); CHKERRQ(ierr);}
	if(temp!=PETSC_NULL){ierr = VecDestroy(&temp); CHKERRQ(ierr);}
	if(lambda!=PETSC_NULL) {ierr = VecDestroy(&lambda); CHKERRQ(ierr);}
	if(r2!=PETSC_NULL)  {ierr = VecDestroy(&r2); CHKERRQ(ierr);}
	if(rhs2!=PETSC_NULL){ierr = VecDestroy(&rhs2); CHKERRQ(ierr);}

	if(uMapping!=PETSC_NULL){ierr = VecDestroy(&uMapping); CHKERRQ(ierr);}
	if(vMapping!=PETSC_NULL){ierr = VecDestroy(&vMapping); CHKERRQ(ierr);}
	if(wMapping!=PETSC_NULL){ierr = VecDestroy(&wMapping); CHKERRQ(ierr);}
	if(pMapping!=PETSC_NULL){ierr = VecDestroy(&pMapping); CHKERRQ(ierr);}

	if(MHat!=PETSC_NULL){ierr = VecDestroy(&MHat); CHKERRQ(ierr);}
	if(RInv!=PETSC_NULL){ierr = VecDestroy(&RInv); CHKERRQ(ierr);}
	if(BN!=PETSC_NULL)  {ierr = VecDestroy(&BN); CHKERRQ(ierr);}

	// Mats
	if(A!=PETSC_NULL)    {ierr = MatDestroy(&A); CHKERRQ(ierr);}
	if(QT!=PETSC_NULL)   {ierr = MatDestroy(&QT); CHKERRQ(ierr);}
	if(BNQ!=PETSC_NULL)  {ierr = MatDestroy(&BNQ); CHKERRQ(ierr);}
	if(QTBNQ!=PETSC_NULL){ierr = MatDestroy(&QTBNQ); CHKERRQ(ierr);}

	// KSPs
	if(ksp1!=PETSC_NULL){ierr = KSPDestroy(&ksp1); CHKERRQ(ierr);}
	if(ksp2!=PETSC_NULL){ierr = KSPDestroy(&ksp2); CHKERRQ(ierr);}

	// Print performance summary to file
	PetscViewer viewer;
	std::string performanceSummaryFileName = caseFolder + "/performanceSummary.txt";
	ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, performanceSummaryFileName.c_str(), &viewer); CHKERRQ(ierr);
	ierr = PetscLogView(viewer); CHKERRQ(ierr);

	return 0;
}

template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::generateRHS1()
{
	PetscErrorCode ierr;
	ierr = VecWAXPY(rhs1, 1.0, rn, bc1); CHKERRQ(ierr);
	ierr = VecPointwiseMult(rhs1, MHat, rhs1); CHKERRQ(ierr);

	return 0;
}

template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::generateRHS2()
{
	PetscErrorCode ierr;
	ierr = VecScale(r2, -1.0); CHKERRQ(ierr);
	ierr = MatMultAdd(QT, qStar, r2, rhs2); CHKERRQ(ierr);

	return 0;
}

template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::stepTime()
{
	PetscErrorCode ierr;

	ierr = PetscLogStagePush(stageSolveIntermediateVelocity); CHKERRQ(ierr);
	ierr = calculateExplicitTerms(); CHKERRQ(ierr);
	ierr = updateBoundaryGhosts(); CHKERRQ(ierr);
	ierr = generateBC1(); CHKERRQ(ierr);
	ierr = generateRHS1(); CHKERRQ(ierr);
	ierr = solveIntermediateVelocity(); CHKERRQ(ierr);
	ierr = PetscLogStagePop(); CHKERRQ(ierr);

	ierr = PetscLogStagePush(stageSolvePoissonSystem); CHKERRQ(ierr);
	ierr = generateR2(); CHKERRQ(ierr);
	ierr = generateRHS2(); CHKERRQ(ierr);
	ierr = solvePoissonSystem(); CHKERRQ(ierr);
	ierr = PetscLogStagePop(); CHKERRQ(ierr);

	ierr = PetscLogStagePush(stageProjectionStep); CHKERRQ(ierr);
	ierr = projectionStep(); CHKERRQ(ierr);
	ierr = PetscLogStagePop(); CHKERRQ(ierr);
	timeStep++;

	return 0;
}

template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::solveIntermediateVelocity()
{
	PetscErrorCode ierr;
	ierr = KSPSolve(ksp1, rhs1, qStar); CHKERRQ(ierr);

	return 0;
}

template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::solvePoissonSystem()
{
	PetscErrorCode ierr;
	ierr = KSPSolve(ksp2, rhs2, lambda); CHKERRQ(ierr);

	return 0;
}

template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::projectionStep()
{
	PetscErrorCode ierr;
	ierr = MatMult(BNQ, lambda, temp); CHKERRQ(ierr);
	ierr = VecWAXPY(q, -1.0, temp, qStar); CHKERRQ(ierr);

	return 0;
}

template <PetscInt dim>
PetscBool NavierStokesSolver<dim>::savePoint()
{
	return (timeStep%simParams->nsave == 0)? PETSC_TRUE : PETSC_FALSE;
}

template <PetscInt dim>
PetscBool NavierStokesSolver<dim>::finished()
{
	return (timeStep >= simParams->nt)? PETSC_TRUE : PETSC_FALSE;
}

template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::generateQTBNQ()
{
	PetscErrorCode ierr;
	PetscInt       pStart;
	PetscLogEvent  GENERATE_QTBNQ;
	
	ierr = PetscLogEventRegister("generateQTBNQ", 0, &GENERATE_QTBNQ); CHKERRQ(ierr);
	ierr = PetscLogEventBegin(GENERATE_QTBNQ, 0, 0, 0, 0); CHKERRQ(ierr);

	ierr = MatMatMult(QT, BNQ, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &QTBNQ); CHKERRQ(ierr);
	
	ierr = VecGetOwnershipRange(lambda, &pStart, NULL); CHKERRQ(ierr);
	if(pStart==0)
	{
		ierr = MatSetValue(QTBNQ, 0, 0, 1.0, ADD_VALUES); CHKERRQ(ierr);
	}
	ierr = MatAssemblyBegin(QTBNQ, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(QTBNQ, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	
	ierr = PetscPrintf(PETSC_COMM_WORLD, "Generated QTBNQ!\n");
	
	ierr = PetscLogEventEnd(GENERATE_QTBNQ, 0, 0, 0, 0); CHKERRQ(ierr);

	return 0;
}

template <PetscInt dim>
void NavierStokesSolver<dim>::countNumNonZeros(PetscInt *cols, size_t numCols, PetscInt rowStart, PetscInt rowEnd, PetscInt &d_nnz, PetscInt &o_nnz)
{
	d_nnz = 0;
	o_nnz = 0;
	for(size_t i=0; i<numCols; i++)
	{
		(cols[i]>=rowStart && cols[i]<rowEnd)? d_nnz++ : o_nnz++;
	}
}

#include "NavierStokes/createDMs.inl"
#include "NavierStokes/createVecs.inl"
#include "NavierStokes/createKSPs.inl"
#include "NavierStokes/createLocalToGlobalMappingsFluxes.inl"
#include "NavierStokes/createLocalToGlobalMappingsLambda.inl"
#include "NavierStokes/initializeMeshSpacings.inl"
#include "NavierStokes/initializeFluxes.inl"
#include "NavierStokes/readFluxes.inl"
#include "NavierStokes/initializeLambda.inl"
#include "NavierStokes/updateBoundaryGhosts.inl"
#include "NavierStokes/calculateExplicitTerms.inl"
#include "NavierStokes/generateDiagonalMatrices.inl"
#include "NavierStokes/generateA.inl"
#include "NavierStokes/generateBC1.inl"
#include "NavierStokes/generateBNQ.inl"
#include "NavierStokes/generateR2.inl"
#include "NavierStokes/writeSimulationInfo.inl"
#include "NavierStokes/writeGrid.inl"
#include "NavierStokes/writeFluxes.inl"
#include "NavierStokes/writeLambda.inl"
#include "NavierStokes/writeData.inl"

template class NavierStokesSolver<2>;
template class NavierStokesSolver<3>;
