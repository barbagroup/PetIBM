/***************************************************************************//**
 * \file diffusiveTermTest.cpp
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \brief Unit-test for the diffusive term.
 */


#include "DiffusiveTerm.h"
#include <navierStokes/NavierStokesSolver.h>
#include <CartesianMesh.h>
#include <FlowDescription.h>
#include <SimulationParameters.h>

#include <memory>

#include <petscksp.h>

#ifndef DIMENSIONS
#define DIMENSIONS 2
#endif


int main(int argc, char **argv)
{
  PetscErrorCode ierr;
  const PetscInt dim = DIMENSIONS;

  ierr = PetscInitialize(&argc, &argv, NULL, NULL); CHKERRQ(ierr);

  char caseFolder[PETSC_MAX_PATH_LEN];
  ierr = PetscOptionsGetString(NULL, "-caseFolder", caseFolder, sizeof(caseFolder), NULL); CHKERRQ(ierr);

  std::string folder(caseFolder);
  FlowDescription FD(folder+"/flowDescription.yaml");
  CartesianMesh CM(folder+"/cartesianMesh.yaml");
  SimulationParameters SP(folder+"/simulationParameters.yaml");

  std::unique_ptr< DiffusiveTerm<dim> > solver(new DiffusiveTerm<dim>(folder, &FD, &SP, &CM));

  ierr = solver->initialize(); CHKERRQ(ierr);

  ierr = solver-> calculateExplicitTerms(); CHKERRQ(ierr);

  ierr = solver->calculateExactSolution(); CHKERRQ(ierr);
  ierr = solver->calculateRelativeError(); CHKERRQ(ierr);
  ierr = solver->writeRelativeError(); CHKERRQ(ierr);

  ierr = solver->finalize(); CHKERRQ(ierr);

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return ierr;
}