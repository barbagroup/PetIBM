#include <string>
#include "petscsys.h"
#include "classes/CartesianMesh.h"
#include "classes/SimulationParameters.h"

int main(int argc,char **argv)
{
	PetscErrorCode       ierr;
	char                 caseFolder[PETSC_MAX_PATH_LEN];
	const PetscInt       dim = 2;
	
	ierr = PetscInitialize(&argc, &argv, NULL, NULL); CHKERRQ(ierr);
	
	ierr = PetscOptionsGetString(NULL, "-caseFolder", caseFolder, sizeof(caseFolder), NULL); CHKERRQ(ierr);
	
	std::string          folder(caseFolder);
	CartesianMesh<dim>   CM("./cases/" + folder + "/cartesianMesh.yaml");
	SimulationParameters SP("./cases/" + folder + "/simulationParameters.yaml");
	
	ierr = PetscFinalize(); CHKERRQ(ierr);
	return 0;
}
