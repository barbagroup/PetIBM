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
  PetscErrorCode ierr;
  const PetscInt dim = DIMENSIONS;
  char           caseFolder[PETSC_MAX_PATH_LEN];
  
  ierr = PetscInitialize(&argc, &argv, NULL, NULL); CHKERRQ(ierr);

  ierr = PetscOptionsGetString(NULL, "-caseFolder", caseFolder, sizeof(caseFolder), NULL); CHKERRQ(ierr);

  std::string          folder(caseFolder);
  FlowDescription      FD(folder+"/flowDescription.yaml");
  CartesianMesh        CM(folder+"/cartesianMesh.yaml");
  SimulationParameters SP(folder+"/simulationParameters.yaml");

  std::unique_ptr< NavierStokesSolver<dim> > solver = createSolver<dim>(folder, &FD, &SP, &CM);
  
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