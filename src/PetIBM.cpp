/***************************************************************************//**
 * \mainpage PetIBM
 *
 *                A PETSc-based Immersed Boundary Method code
 *
 * \author Anush Krishnan (anush@bu.edu)
 *         Olivier Mesnard (mesnardo@gwu.edu)
 */

/***************************************************************************//**
 * \file PetIBM.cpp
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Main source-file of `PetIBM`.
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

  ierr = PetscPrintf(PETSC_COMM_WORLD, "\n======================\n*** PetIBM - Start ***\n======================\n"); CHKERRQ(ierr);


  // parse command-line to get simulation directory
  char dir[PETSC_MAX_PATH_LEN];
  PetscBool found;
  ierr = PetscOptionsGetString(NULL, NULL, "-directory", dir, sizeof(dir), &found); CHKERRQ(ierr);
  std::string directory(".");
  if (found)
    directory = dir;
  ierr = PetscPrintf(PETSC_COMM_WORLD, "directory: %s\n", directory.c_str()); CHKERRQ(ierr);

  // read different input files
  CartesianMesh mesh(directory+"/cartesianMesh.yaml");
  ierr = mesh.printInfo(); CHKERRQ(ierr);
  ierr = mesh.write(directory+"/grid.txt");
  FlowDescription<dim> flow(directory+"/flowDescription.yaml");
  ierr = flow.printInfo(); CHKERRQ(ierr);
  SimulationParameters parameters(directory, directory+"/simulationParameters.yaml");
  ierr = parameters.printInfo(); CHKERRQ(ierr);

  std::unique_ptr< NavierStokesSolver<dim> > solver = createSolver<dim>(&mesh,
                                                                        &flow,
                                                                        &parameters);
  
  ierr = solver->initialize(); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "\nInitialization complete!\n\n"); CHKERRQ(ierr);
  
  while(!solver->finished())
  {
    ierr = solver->stepTime(); CHKERRQ(ierr);
    ierr = solver->writeData(); CHKERRQ(ierr);
  }
  
  ierr = solver->finalize(); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD, "\n=====================\n*** PetIBM - Done ***\n=====================\n"); CHKERRQ(ierr);

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;
} // main