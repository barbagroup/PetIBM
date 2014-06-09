#include "FlowDescription.h"
#include "CartesianMesh.h"
#include "SimulationParameters.h"
#include "createSolver.h"
#include <petscsys.h>
#include <iostream>
#include <string>
#include <memory>

int main(int argc,char **argv)
{
	PetscErrorCode ierr;
	const PetscInt dim = 2;
	char           caseFolder[PETSC_MAX_PATH_LEN];
	
	ierr = PetscInitialize(&argc, &argv, NULL, NULL); CHKERRQ(ierr);
	
	ierr = PetscOptionsGetString(NULL, "-caseFolder", caseFolder, sizeof(caseFolder), NULL); CHKERRQ(ierr);

	std::string             folder(caseFolder);
	FlowDescription         FD(folder+"/flowDescription.yaml");
	CartesianMesh           CM(folder+"/cartesianMesh.yaml");
	SimulationParameters    SP(folder+"/simulationParameters.yaml");

	std::unique_ptr< NavierStokesSolver<dim> > solver = createSolver<dim>(folder, &FD, &SP, &CM);
	
	ierr = solver->initialise(); CHKERRQ(ierr);
	solver->writeSimulationInfo(folder);
	solver->writeGrid(folder);
	
	while(!solver->finished())
	{
		ierr = solver->stepTime(); CHKERRQ(ierr);
		ierr = solver->writeData(folder); CHKERRQ(ierr);
	}
	ierr = solver->finalise(); CHKERRQ(ierr);

	ierr = PetscFinalize(); CHKERRQ(ierr);
	return 0;
}
