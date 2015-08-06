/***************************************************************************//**
 * \file TairaColoniusSolver.cpp
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the methods of the class \c TairaColoniusSolver.
 */


#include "TairaColoniusSolver.h"

#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <sys/stat.h>

#include "yaml-cpp/yaml.h"
#include <petscdmcomposite.h>


template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::initialize()
{
  PetscErrorCode ierr;

  ierr = PetscLogStagePush(NavierStokesSolver<dim>::stageInitialize); CHKERRQ(ierr);
  ierr = initializeBodies(); CHKERRQ(ierr);
  ierr = calculateCellIndices(); CHKERRQ(ierr);
  ierr = createDMs(); CHKERRQ(ierr);
  ierr = createGlobalMappingBodies(); CHKERRQ(ierr);
  ierr = NavierStokesSolver<dim>::initializeCommon(); CHKERRQ(ierr);
  ierr = PetscLogStagePop(); CHKERRQ(ierr);

  return 0;
} // initialize

template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::finalize()
{
  PetscErrorCode ierr;

  ierr = NavierStokesSolver<dim>::finalize();

  // DMs
  if(bda!=PETSC_NULL) {ierr = DMDestroy(&bda); CHKERRQ(ierr);}
  // Mats
  if(ET!=PETSC_NULL)  {ierr = MatDestroy(&ET); CHKERRQ(ierr);}
  // Vecs
  if(regularizedForce!=PETSC_NULL){ierr = VecDestroy(&regularizedForce); CHKERRQ(ierr);}
  if(nullSpaceVec!=PETSC_NULL){ierr = VecDestroy(&nullSpaceVec); CHKERRQ(ierr);}

  return 0;
}  // finalize

template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::createDMs()
{
  PetscErrorCode ierr;
  ierr = NavierStokesSolver<dim>::createDMs(); CHKERRQ(ierr); 
  ierr = generateBodyInfo(); CHKERRQ(ierr);
  ierr = DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, x.size(), dim, 0, &numBoundaryPointsOnProcess.front(), &bda); CHKERRQ(ierr);
  ierr = DMCompositeAddDM(NavierStokesSolver<dim>::lambdaPack, bda); CHKERRQ(ierr);

  return 0;
} // createDMs

template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::createVecs()
{
  PetscErrorCode ierr;

  ierr = NavierStokesSolver<dim>::createVecs();
  ierr = VecDuplicate(NavierStokesSolver<dim>::q, &regularizedForce); CHKERRQ(ierr);
  ierr = VecDuplicate(NavierStokesSolver<dim>::lambda, &nullSpaceVec); CHKERRQ(ierr);

  return 0;
} // createVecs

template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::writeData()
{
  NavierStokesSolver<dim>::writeData();
  calculateForce();
  writeForces();

  return 0;
} // writeData

template <PetscInt dim>
PetscReal TairaColoniusSolver<dim>::dhRoma(PetscReal x, PetscReal h)
{
  PetscReal r = fabs(x)/h;
  if(r>1.5) return 0.0;
  if(r>0.5 && r<=1.5) return 1.0/(6*h)*( 5.0 - 3.0*r - sqrt(-3.0*(1-r)*(1-r) + 1.0) );
  return 1.0/(3*h)*( 1.0 + sqrt(-3.0*r*r + 1.0) );
} // dhRoma

template <PetscInt dim>
PetscReal TairaColoniusSolver<dim>::delta(PetscReal x, PetscReal y, PetscReal h)
{
  return dhRoma(x, h) * dhRoma(y, h);
} // delta

template <PetscInt dim>
PetscReal TairaColoniusSolver<dim>::delta(PetscReal x, PetscReal y, PetscReal z, PetscReal h)
{
  return dhRoma(x, h) * dhRoma(y, h) * dhRoma(z, h);
} // delta

#include "inline/setNullSpace.inl"
#include "inline/calculateCellIndices.inl"
#include "inline/generateBodyInfo.inl"
#include "inline/generateBNQ.inl"
#include "inline/generateR2.inl"
#include "inline/initializeBodies.inl"
#include "inline/createGlobalMappingBodies.inl"
#include "inline/isInfluenced.inl"
#include "inline/calculateForce.inl"
#include "inline/io.inl"

template class TairaColoniusSolver<2>;
template class TairaColoniusSolver<3>;