#include "createSolver.h"
#include "gtest/gtest.h"
#include <petscdmcomposite.h>

class NavierStokesTest : public ::testing::Test
{
public:
	std::string           folder;
	FlowDescription       FD;
	CartesianMesh         CM;
	SimulationParameters  SP;
	std::unique_ptr< NavierStokesSolver<2> > solver;
	Vec                   lambdaGold, error;

	NavierStokesTest()
	{
		char           caseFolder[PETSC_MAX_PATH_LEN];

		lambdaGold = PETSC_NULL;
		
		// read case folder
		PetscOptionsGetString(NULL, "-caseFolder", caseFolder, sizeof(caseFolder), NULL);
		
		// read input files and create solver
		folder = std::string(caseFolder);
		FD = FlowDescription(folder+"/flowDescription.yaml");
		CM = CartesianMesh(folder+"/cartesianMesh.yaml");
		SP = SimulationParameters(folder+"/simulationParameters.yaml");
		solver = createSolver<2>(folder, &FD, &SP, &CM);
	}

	virtual void SetUp()
	{
		// initialise solver
		solver->initialize();

		// perform the simulation
		while(!solver->finished())
		{
			solver->stepTime();
		}
	}

	virtual void TearDown()
	{
		solver->finalize();
		if(lambdaGold!=PETSC_NULL) VecDestroy(&lambdaGold);
	}
};

TEST_F(NavierStokesTest, ComparePhi)
{
	PetscViewer viewer;
	Vec         phi;
	PetscReal   errorNorm, goldNorm;

	// create vector to store the gold data and the error
	VecDuplicate(solver->lambda, &lambdaGold);
	VecDuplicate(solver->lambda, &error);
	DMCreateGlobalVector(solver->lambdaPack, &lambdaGold);

	// read the gold data from file
	DMCompositeGetAccess(solver->lambdaPack, lambdaGold, &phi);
	std::string fileName = folder + "/phi.dat";
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, fileName.c_str(), FILE_MODE_READ, &viewer);
	VecLoad(phi, viewer);
	PetscViewerDestroy(&viewer);
	DMCompositeRestoreAccess(solver->lambdaPack, lambdaGold, &phi);

	// print the difference
	VecWAXPY(error, -1, solver->lambda, lambdaGold);
	VecNorm(error, NORM_2, &errorNorm);
	VecNorm(lambdaGold, NORM_2, &goldNorm);
	//PetscPrintf(PETSC_COMM_WORLD, "%f\n", errorNorm/goldNorm);

	EXPECT_LT(errorNorm/goldNorm, 5e-4);
}

int main(int argc, char **argv)
{
	PetscErrorCode ierr, result;

	::testing::InitGoogleTest(&argc, argv);
	ierr = PetscInitialize(&argc, &argv, NULL, NULL); CHKERRQ(ierr);
	result = RUN_ALL_TESTS();
	ierr = PetscFinalize(); CHKERRQ(ierr);

	return result;
}
