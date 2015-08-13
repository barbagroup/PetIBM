/***************************************************************************//**
 * \mainpage PetIBM
 *
 *                A PETSc-based Immersed Boundary Method code
 *
 * \author Anush Krishnan (anush@bu.edu)
 */

/***************************************************************************//**
 * \file PetIBM.cpp
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Main source-file of \c PetIBM.
 */


#include "createSolver.h"

#ifndef DIMENSIONS
#define DIMENSIONS 2
#endif


int main(int argc,char **argv)
{
  const PetscInt dim = DIMENSIONS;
  PetscErrorCode ierr;
  
  ierr = PetscInitialize(&argc, &argv, NULL, NULL); CHKERRQ(ierr);

  // parse command-line to get simulation directory
  char dir[PETSC_MAX_PATH_LEN];
  ierr = PetscOptionsGetString(NULL, "-directory", dir, sizeof(dir), NULL); CHKERRQ(ierr);
  std::string directory(dir);

  // read different input files
  CartesianMesh cartesianMesh(directory);
  FlowDescription<dim> flowDescription(directory);
  SimulationParameters simulationParameters(directory);

  std::unique_ptr< NavierStokesSolver<dim> > solver = createSolver<dim>(&cartesianMesh,
                                                                        &flowDescription, 
                                                                        &simulationParameters);
  
  ierr = solver->initialize(); CHKERRQ(ierr);
  
  while(!solver->finished())
  {
    ierr = solver->stepTime(); CHKERRQ(ierr);
    ierr = solver->writeData(); CHKERRQ(ierr);
  }
  
  ierr = solver->finalize(); CHKERRQ(ierr);

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;
} // main