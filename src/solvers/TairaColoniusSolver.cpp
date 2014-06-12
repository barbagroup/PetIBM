#include "TairaColoniusSolver.h"
#include <petscsys.h>
#include <petscdmcomposite.h>
#include <iostream>
#include <string>

#include "TairaColonius/createDMs.inl"
#include "TairaColonius/generateBNQ.inl"
#include "TairaColonius/generateET.inl"
#include "TairaColonius/generateR2.inl"
#include "TairaColonius/initialiseBodies.inl"
#include "TairaColonius/createGlobalMappingBodies.inl"
#include "TairaColonius/calculateForce.inl"
#include "TairaColonius/writeForces.inl"

template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::initialise()
{
	PetscErrorCode ierr;

	initialiseBodies();
	ierr = createDMs(); CHKERRQ(ierr);
	ierr = NavierStokesSolver<dim>::createVecs(); CHKERRQ(ierr);
	ierr = VecDuplicate(NavierStokesSolver<dim>::q, &temp); CHKERRQ(ierr);

	NavierStokesSolver<dim>::initialiseMeshSpacings();
	ierr = NavierStokesSolver<dim>::initialiseFluxes(); CHKERRQ(ierr);
	ierr = NavierStokesSolver<dim>::updateBoundaryGhosts(); CHKERRQ(ierr);

	ierr = createGlobalMappingBodies(); CHKERRQ(ierr);
	ierr = NavierStokesSolver<dim>::createLocalToGlobalMappingsFluxes(); CHKERRQ(ierr);
	ierr = NavierStokesSolver<dim>::createLocalToGlobalMappingsLambda(); CHKERRQ(ierr);

	ierr = NavierStokesSolver<dim>::generateDiagonalMatrices(); CHKERRQ(ierr);
	ierr = NavierStokesSolver<dim>::generateA(); CHKERRQ(ierr);
	ierr = generateBNQ(); CHKERRQ(ierr);
	ierr = NavierStokesSolver<dim>::generateQTBNQ(); CHKERRQ(ierr);
	ierr = NavierStokesSolver<dim>::createKSPs(); CHKERRQ(ierr);
	ierr = generateET(); CHKERRQ(ierr);

	return 0;
}

template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::finalise()
{
	PetscErrorCode ierr;

	ierr = NavierStokesSolver<dim>::finalise();

	// DMs
	if(bda!=PETSC_NULL) {ierr = DMDestroy(&bda); CHKERRQ(ierr);}
	// Mats
	if(ET!=PETSC_NULL)  {ierr = MatDestroy(&ET); CHKERRQ(ierr);}
	// Vecs
	if(temp!=PETSC_NULL){ierr = VecDestroy(&temp); CHKERRQ(ierr);}

	return 0;
}

template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::writeData()
{
	NavierStokesSolver<dim>::writeData();
	calculateForce();
	writeForces();

	return 0;
}

template <PetscInt dim>
PetscReal TairaColoniusSolver<dim>::dhRoma(PetscReal x, PetscReal h)
{
	PetscReal r = fabs(x)/h;
	if(r>1.5)
		return 0.0;
	else if(r>0.5 && r<=1.5)
		return 1.0/(6*h)*( 5.0 - 3.0*r - sqrt(-3.0*(1-r)*(1-r) + 1.0) );
	else
		return 1.0/(3*h)*( 1.0 + sqrt(-3.0*r*r + 1.0) );
}

template <PetscInt dim>
PetscReal TairaColoniusSolver<dim>::delta(PetscReal x, PetscReal y, PetscReal h)
{
	return dhRoma(x, h) * dhRoma(y, h);
}

template <PetscInt dim>
PetscReal TairaColoniusSolver<dim>::delta(PetscReal x, PetscReal y, PetscReal z, PetscReal h)
{
	return dhRoma(x, h) * dhRoma(y, h) * dhRoma(z, h);
}

template <PetscInt dim>
bool TairaColoniusSolver<dim>::isInfluenced(PetscReal xGrid, PetscReal yGrid, PetscReal xBody, PetscReal yBody, PetscReal radius, PetscReal *disp)
{
	PetscReal width[2];
	PetscReal nx = NavierStokesSolver<dim>::mesh->nx,
	          ny = NavierStokesSolver<dim>::mesh->ny;
	
	std::vector<PetscReal> &xMesh = NavierStokesSolver<dim>::mesh->x,
	                       &yMesh = NavierStokesSolver<dim>::mesh->y;

	width[0] = xMesh[nx]-xMesh[0];
	width[1] = yMesh[ny]-yMesh[0];

	disp[0] = fabs(xGrid - xBody);
	disp[1] = fabs(yGrid - yBody);

	if(NavierStokesSolver<dim>::flowDesc->bc[0][XPLUS].type==PERIODIC && disp[0]>width[0]-disp[0]) disp[0] = width[0] - disp[0];
	if(NavierStokesSolver<dim>::flowDesc->bc[0][YPLUS].type==PERIODIC && disp[1]>width[1]-disp[1]) disp[1] = width[1] - disp[1];

	return (disp[0] < radius && disp[1] < radius);
}

template <PetscInt dim>
bool TairaColoniusSolver<dim>::isInfluenced(PetscReal xGrid, PetscReal yGrid, PetscReal zGrid, PetscReal xBody, PetscReal yBody, PetscReal zBody, PetscReal radius, PetscReal *disp)
{
	PetscReal width[3];
	PetscReal nx = NavierStokesSolver<dim>::mesh->nx,
	          ny = NavierStokesSolver<dim>::mesh->ny,
	          nz = NavierStokesSolver<dim>::mesh->nz;

	std::vector<PetscReal> &xMesh = NavierStokesSolver<dim>::mesh->x,
	                       &yMesh = NavierStokesSolver<dim>::mesh->y,
	                       &zMesh = NavierStokesSolver<dim>::mesh->z;

	width[0] = xMesh[nx]-xMesh[0];
	width[1] = yMesh[ny]-yMesh[0];
	width[2] = zMesh[nz]-zMesh[0];

	disp[0] = fabs(xGrid - xBody);
	disp[1] = fabs(yGrid - yBody);
	disp[2] = fabs(zGrid - zBody);

	if(NavierStokesSolver<dim>::flowDesc->bc[0][XPLUS].type==PERIODIC && disp[0]>width[0]-disp[0]) disp[0] = width[0] - disp[0];
	if(NavierStokesSolver<dim>::flowDesc->bc[0][YPLUS].type==PERIODIC && disp[1]>width[1]-disp[1]) disp[1] = width[1] - disp[1];
	if(NavierStokesSolver<dim>::flowDesc->bc[0][ZPLUS].type==PERIODIC && disp[2]>width[2]-disp[2]) disp[2] = width[2] - disp[2];

	return (disp[0] < radius && disp[1] < radius && disp[2] < radius);
}

template class TairaColoniusSolver<2>;
template class TairaColoniusSolver<3>;