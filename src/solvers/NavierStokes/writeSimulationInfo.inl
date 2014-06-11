#include <fstream>

template <>
void NavierStokesSolver<2>::writeSimulationInfo()
{
	PetscInt    rank;
	
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
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
}

template <>
void NavierStokesSolver<3>::writeSimulationInfo()
{
	PetscInt    rank;
	
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
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
}