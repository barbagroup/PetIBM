#include "FlowDescription.h"
#include "CartesianMesh.h"
#include "SimulationParameters.h"
#include "NavierStokesSolver.h"
#include <petscsys.h>
#include <string>

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
	NavierStokesSolver<dim> *solver = NULL;
	
	solver = NavierStokesSolver<dim>::createSolver(FD, SP, CM);
	
	//solver->initialise();
	
	if(!solver)
		delete solver;
	ierr = PetscFinalize(); CHKERRQ(ierr);
	return 0;
}
