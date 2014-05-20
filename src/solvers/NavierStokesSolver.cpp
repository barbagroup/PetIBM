#include "NavierStokesSolver.h"
#include <petscsys.h>
#include <iostream>
#include <string>

template <PetscInt dim>
void NavierStokesSolver<dim>::initialise()
{
	createDMs();
	createVecs();
	createLocalToGlobalMappings();
	initialiseMeshSpacings();
	initialiseFluxes();
	updateBoundaryGhosts();

	generateDiagonalMatrices();
	generateA();
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
	if(pack!=PETSC_NULL){ierr = DMDestroy(&pack); CHKERRV(ierr);}
	
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

	if(uMapping!=PETSC_NULL){ierr = VecDestroy(&uMapping); CHKERRV(ierr);}
	if(vMapping!=PETSC_NULL){ierr = VecDestroy(&vMapping); CHKERRV(ierr);}
	if(wMapping!=PETSC_NULL){ierr = VecDestroy(&wMapping); CHKERRV(ierr);}

	if(MHat!=PETSC_NULL){ierr = VecDestroy(&MHat); CHKERRV(ierr);}
	if(RInv!=PETSC_NULL){ierr = VecDestroy(&RInv); CHKERRV(ierr);}
	if(BN!=PETSC_NULL)  {ierr = VecDestroy(&BN); CHKERRV(ierr);}

	// Mats
	if(A!=PETSC_NULL)    {ierr = MatDestroy(&A); CHKERRV(ierr);}
	if(QT!=PETSC_NULL)   {ierr = MatDestroy(&QT); CHKERRV(ierr);}
	if(BNQ!=PETSC_NULL)  {ierr = MatDestroy(&BNQ); CHKERRV(ierr);}
	if(QTBNQ!=PETSC_NULL){ierr = MatDestroy(&QTBNQ); CHKERRV(ierr);}
}

template <PetscInt dim>
void NavierStokesSolver<dim>::generateRHS1()
{
	PetscErrorCode ierr;
	ierr = VecWAXPY(rhs1, 1.0, rn, bc1); CHKERRV(ierr);
}

template <PetscInt dim>
void NavierStokesSolver<dim>::stepTime()
{
	timeStep++;
}

template <PetscInt dim>
void NavierStokesSolver<dim>::writeData()
{
}

template <PetscInt dim>
bool NavierStokesSolver<dim>::finished()
{
	return (timeStep < simParams->nt)? false : true;
}

#include "NavierStokes/createDMs.inl"
#include "NavierStokes/createVecs.inl"
#include "NavierStokes/createLocalToGlobalMappings.inl"
#include "NavierStokes/initialiseMeshSpacings.inl"
#include "NavierStokes/initialiseFluxes.inl"
#include "NavierStokes/updateBoundaryGhosts.inl"
#include "NavierStokes/calculateExplicitTerms.inl"
#include "NavierStokes/generateDiagonalMatrices.inl"
#include "NavierStokes/generateA.inl"
#include "NavierStokes/generateBC1.inl"

template class NavierStokesSolver<2>;
template class NavierStokesSolver<3>;