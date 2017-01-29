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
  ierr = createDMs(); CHKERRQ(ierr);
  ierr = createGlobalMappingBodies(); CHKERRQ(ierr);
  ierr = NavierStokesSolver<dim>::initializeCommon(); CHKERRQ(ierr);

  ierr = PetscLogStagePop(); CHKERRQ(ierr);

  return 0;
} // initialize


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
#include "inline/calculateForces.inl"
#include "inline/io.inl"


// dimensions specialization
template class TairaColoniusSolver<2>;
template class TairaColoniusSolver<3>;