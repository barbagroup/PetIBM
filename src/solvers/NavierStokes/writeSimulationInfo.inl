#include <fstream>

template <>
PetscErrorCode NavierStokesSolver<2>::writeSimulationInfo()
{
	PetscErrorCode ierr;
	PetscInt       rank;
	
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
	
	if(rank==0)
	{
		std::ofstream f(caseFolder+"/simulationInfo.txt");
		f << "-nx\t" << mesh->nx << '\n';
		f << "-ny\t" << mesh->ny << '\n';
		f << "-nt\t" << simParams->nt << '\n';
		f << "-nsave\t" << simParams->nsave << '\n';
		(flowDesc->bc[0][XPLUS].type==PERIODIC)? f << "-xperiodic\tTrue\n" : f << "-xperiodic\tFalse\n";
		(flowDesc->bc[0][YPLUS].type==PERIODIC)? f << "-yperiodic\tTrue\n" : f << "-yperiodic\tFalse\n";
		f.close();
	}

	return 0;
}

template <>
PetscErrorCode NavierStokesSolver<3>::writeSimulationInfo()
{
	PetscErrorCode ierr;
	PetscInt       rank;
	
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
	
	if(rank==0)
	{
		std::ofstream f(caseFolder+"/simulationInfo.txt");
		f << "-nx\t" << mesh->nx << '\n';
		f << "-ny\t" << mesh->ny << '\n';
		f << "-nz\t" << mesh->nz << '\n';
		f << "-nt\t" << simParams->nt << '\n';
		f << "-nsave\t" << simParams->nsave << '\n';
		(flowDesc->bc[0][XPLUS].type==PERIODIC)? f << "-xperiodic\tTrue\n" : f << "-xperiodic\tFalse\n";
		(flowDesc->bc[0][YPLUS].type==PERIODIC)? f << "-yperiodic\tTrue\n" : f << "-yperiodic\tFalse\n";
		(flowDesc->bc[0][ZPLUS].type==PERIODIC)? f << "-zperiodic\tTrue\n" : f << "-zperiodic\tFalse\n";
		f.close();
	}

	return 0;
}