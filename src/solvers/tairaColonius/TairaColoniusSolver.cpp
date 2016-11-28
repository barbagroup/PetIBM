/***************************************************************************//**
 * \file TairaColoniusSolver.cpp
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the methods of the class `TairaColoniusSolver`.
 */


#include "TairaColoniusSolver.h"

#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <sys/stat.h>

#include <petscdmcomposite.h>


/**
 * \brief Constructor. Calls parent's contructor and initializes additioner pointers.
 */
template <PetscInt dim>
TairaColoniusSolver<dim>::TairaColoniusSolver(CartesianMesh *cartesianMesh, 
                                              FlowDescription<dim> *flowDescription, 
                                              SimulationParameters *simulationParameters)
                        : NavierStokesSolver<dim>::NavierStokesSolver(cartesianMesh, 
                                                                      flowDescription, 
                                                                      simulationParameters)
{
  bda = PETSC_NULL;
  ET = PETSC_NULL;
  nullSpaceVec = PETSC_NULL;
  regularizedForce = PETSC_NULL;
} // TairaColoniusSolver


/**
 * \brief Initializes the solver.
 */
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


/**
 * \brief Gets the indices of cells owning a body point.
 */
template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::calculateCellIndices()
{
  PetscErrorCode ierr;

  for (PetscInt l=0; l<numBodies; l++)
  {
    ierr = bodies[l].getCellOwners(NavierStokesSolver<dim>::mesh); CHKERRQ(ierr);
  }

  return 0;
} // calculateCellIndices


/**
 * \brief Discrete delta function from Roma et al. (1999).
 */
template <PetscInt dim>
PetscReal TairaColoniusSolver<dim>::dhRoma(PetscReal x, PetscReal h)
{
  PetscReal r = fabs(x)/h;
  if (r > 1.5)
    return 0.0;
  if (r > 0.5 && r <= 1.5)
    return 1.0/(6*h)*( 5.0 - 3.0*r - sqrt(-3.0*(1-r)*(1-r) + 1.0) );
  return 1.0/(3*h)*( 1.0 + sqrt(-3.0*r*r + 1.0) );
} // dhRoma


/**
 * \brief Two-dimensional discrete delta function.
 */
template <PetscInt dim>
PetscReal TairaColoniusSolver<dim>::delta(PetscReal x, PetscReal y, PetscReal h)
{
  return dhRoma(x, h) * dhRoma(y, h);
} // delta


/**
 * \brief Three-dimensional discrete delta function.
 */
template <PetscInt dim>
PetscReal TairaColoniusSolver<dim>::delta(PetscReal x, PetscReal y, PetscReal z, PetscReal h)
{
  return dhRoma(x, h) * dhRoma(y, h) * dhRoma(z, h);
} // delta


/**
 * \brief Destroys PETSc objects.
 */
template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::finalize()
{
  PetscErrorCode ierr;

  ierr = NavierStokesSolver<dim>::finalize();
  // DMs
  if (bda != PETSC_NULL) {ierr = DMDestroy(&bda); CHKERRQ(ierr);}
  // Mats
  if (ET != PETSC_NULL)  {ierr = MatDestroy(&ET); CHKERRQ(ierr);}
  // Vecs
  if (regularizedForce != PETSC_NULL){ierr = VecDestroy(&regularizedForce); CHKERRQ(ierr);}
  if (nullSpaceVec != PETSC_NULL)    {ierr = VecDestroy(&nullSpaceVec); CHKERRQ(ierr);}

  return 0;
}  // finalize


#include "inline/createDMs.inl"
#include "inline/createVecs.inl"
#include "inline/setNullSpace.inl"
#include "inline/generateBodyInfo.inl"
#include "inline/generateBNQ.inl"
#include "inline/generateR2.inl"
#include "inline/initializeBodies.inl"
#include "inline/createGlobalMappingBodies.inl"
#include "inline/isInfluenced.inl"
#include "inline/calculateForce.inl"
#include "inline/io.inl"


// dimensions specialization
template class TairaColoniusSolver<2>;
template class TairaColoniusSolver<3>;