#include <sys/stat.h>
#include <sstream>

template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::writeData()
{
	PetscErrorCode  ierr;
	PetscInt        rank;

	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

	if(rank==0)
	{
		PetscInt its1, its2;
		std::string filename = caseFolder + "/iterationCount.txt";
		if(timeStep==1)
		{
			iterationsFile.open(filename.c_str());
		}
		else
		{
			iterationsFile.open(filename.c_str(), std::ios::out | std::ios::app);
		}
		ierr = KSPGetIterationNumber(ksp1, &its1); CHKERRQ(ierr);
		ierr = KSPGetIterationNumber(ksp2, &its2); CHKERRQ(ierr);
		iterationsFile << timeStep << '\t' << its1 << '\t' << its2 << std::endl;
		iterationsFile.close();
	}

	if(timeStep%simParams->nsave == 0)
	{
		ierr = writeFluxes(); CHKERRQ(ierr);
		ierr = writeLambda(); CHKERRQ(ierr);
	}
	
	return 0;
}