#include "NavierStokesSolver.h"
#include <petscsys.h>
#include <petscdmcomposite.h>
#include <iostream>
#include <string>

template <PetscInt dim>
void NavierStokesSolver<dim>::initialise()
{
	createDMs();
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

template <PetscInt dim>
void NavierStokesSolver<dim>::finalise()
{
	PetscErrorCode ierr;
	
	// DMs
	if(pda!=PETSC_NULL) {ierr = DMDestroy(&pda); CHKERRV(ierr);}
	if(uda!=PETSC_NULL) {ierr = DMDestroy(&uda); CHKERRV(ierr);}
	if(vda!=PETSC_NULL) {ierr = DMDestroy(&vda); CHKERRV(ierr);}
	if(wda!=PETSC_NULL) {ierr = DMDestroy(&wda); CHKERRV(ierr);}
	if(qPack!=PETSC_NULL){ierr = DMDestroy(&qPack); CHKERRV(ierr);}
	if(lambdaPack!=PETSC_NULL){ierr = DMDestroy(&lambdaPack); CHKERRV(ierr);}
	
	// Vecs
	if(q!=PETSC_NULL)    {ierr = VecDestroy(&q); CHKERRV(ierr);}
	if(qStar!=PETSC_NULL){ierr = VecDestroy(&qStar); CHKERRV(ierr);}
	
	if(qxLocal!=PETSC_NULL){ierr = VecDestroy(&qxLocal); CHKERRV(ierr);}
	if(qyLocal!=PETSC_NULL){ierr = VecDestroy(&qyLocal); CHKERRV(ierr);}
	if(qzLocal!=PETSC_NULL){ierr = VecDestroy(&qzLocal); CHKERRV(ierr);}

	if(H!=PETSC_NULL)   {ierr = VecDestroy(&H); CHKERRV(ierr);}
	if(rn!=PETSC_NULL)  {ierr = VecDestroy(&rn); CHKERRV(ierr);}
	if(bc1!=PETSC_NULL) {ierr = VecDestroy(&bc1); CHKERRV(ierr);}
	if(rhs1!=PETSC_NULL){ierr = VecDestroy(&rhs1); CHKERRV(ierr);}
	if(temp!=PETSC_NULL){ierr = VecDestroy(&temp); CHKERRV(ierr);}
	if(lambda!=PETSC_NULL) {ierr = VecDestroy(&lambda); CHKERRV(ierr);}
	if(r2!=PETSC_NULL)  {ierr = VecDestroy(&r2); CHKERRV(ierr);}
	if(rhs2!=PETSC_NULL){ierr = VecDestroy(&rhs2); CHKERRV(ierr);}

	if(uMapping!=PETSC_NULL){ierr = VecDestroy(&uMapping); CHKERRV(ierr);}
	if(vMapping!=PETSC_NULL){ierr = VecDestroy(&vMapping); CHKERRV(ierr);}
	if(wMapping!=PETSC_NULL){ierr = VecDestroy(&wMapping); CHKERRV(ierr);}
	if(pMapping!=PETSC_NULL){ierr = VecDestroy(&pMapping); CHKERRV(ierr);}

	if(MHat!=PETSC_NULL){ierr = VecDestroy(&MHat); CHKERRV(ierr);}
	if(RInv!=PETSC_NULL){ierr = VecDestroy(&RInv); CHKERRV(ierr);}
	if(BN!=PETSC_NULL)  {ierr = VecDestroy(&BN); CHKERRV(ierr);}

	// Mats
	if(A!=PETSC_NULL)    {ierr = MatDestroy(&A); CHKERRV(ierr);}
	if(QT!=PETSC_NULL)   {ierr = MatDestroy(&QT); CHKERRV(ierr);}
	if(BNQ!=PETSC_NULL)  {ierr = MatDestroy(&BNQ); CHKERRV(ierr);}
	if(QTBNQ!=PETSC_NULL){ierr = MatDestroy(&QTBNQ); CHKERRV(ierr);}

	// KSPs
	if(ksp1!=PETSC_NULL){ierr = KSPDestroy(&ksp1); CHKERRV(ierr);}
	if(ksp2!=PETSC_NULL){ierr = KSPDestroy(&ksp2); CHKERRV(ierr);}

	// PCs
	//if(pc2!=PETSC_NULL){ierr = PCDestroy(&pc2); CHKERRV(ierr);}
}

template <PetscInt dim>
void NavierStokesSolver<dim>::generateRHS1()
{
	PetscErrorCode ierr;
	ierr = VecWAXPY(rhs1, 1.0, rn, bc1); CHKERRV(ierr);
	ierr = VecPointwiseMult(rhs1, MHat, rhs1);
}

template <PetscInt dim>
void NavierStokesSolver<dim>::generateRHS2()
{
	PetscErrorCode ierr;
	ierr = VecScale(r2, -1.0); CHKERRV(ierr);
	ierr = MatMultAdd(QT, qStar, r2, rhs2); CHKERRV(ierr);
}

template <PetscInt dim>
void NavierStokesSolver<dim>::stepTime()
{
	calculateExplicitTerms();
	updateBoundaryGhosts();
	generateBC1();
	generateRHS1();
	solveIntermediateVelocity();
	generateR2();
	generateRHS2();
	solvePoissonSystem();
	projectionStep();
	timeStep++;
}

template <PetscInt dim>
void NavierStokesSolver<dim>::solveIntermediateVelocity()
{
	PetscErrorCode ierr;
	ierr = KSPSolve(ksp1, rhs1, qStar); CHKERRV(ierr);
}

template <PetscInt dim>
void NavierStokesSolver<dim>::solvePoissonSystem()
{
	PetscErrorCode ierr;
	ierr = KSPSolve(ksp2, rhs2, lambda); CHKERRV(ierr);
}

template <PetscInt dim>
void NavierStokesSolver<dim>::projectionStep()
{
	PetscErrorCode ierr;
	ierr = MatMult(BNQ, lambda, temp); CHKERRV(ierr);
	ierr = VecWAXPY(q, -1.0, temp, qStar);
}

template <PetscInt dim>
bool NavierStokesSolver<dim>::savePoint()
{
	return (timeStep%simParams->nsave == 0);
}

template <PetscInt dim>
bool NavierStokesSolver<dim>::finished()
{
	return (timeStep < simParams->nt)? false : true;
}

template <PetscInt dim>
void NavierStokesSolver<dim>::generateQTBNQ()
{
	PetscErrorCode ierr;
	PetscInt       pStart;

	ierr = MatMatMult(QT, BNQ, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &QTBNQ); CHKERRV(ierr);
	
	ierr = VecGetOwnershipRange(lambda, &pStart, NULL); CHKERRV(ierr);
	if(pStart==0)
	{
		ierr = MatSetValue(QTBNQ, 0, 0, 1.0, ADD_VALUES); CHKERRV(ierr);
	}
	ierr = MatAssemblyBegin(QTBNQ, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
	ierr = MatAssemblyEnd(QTBNQ, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
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
#include "NavierStokes/initialiseMeshSpacings.inl"
#include "NavierStokes/initialiseFluxes.inl"
#include "NavierStokes/updateBoundaryGhosts.inl"
#include "NavierStokes/calculateExplicitTerms.inl"
#include "NavierStokes/generateDiagonalMatrices.inl"
#include "NavierStokes/generateA.inl"
#include "NavierStokes/generateBC1.inl"
#include "NavierStokes/generateBNQ.inl"
#include "NavierStokes/generateR2.inl"
#include "NavierStokes/writeData.inl"
#include "NavierStokes/writeSimulationInfo.inl"
#include "NavierStokes/writeGrid.inl"

template class NavierStokesSolver<2>;
template class NavierStokesSolver<3>;
