#include "SimulationParameters.h"
#include "yaml-cpp/yaml.h"
#include <fstream>

SolverType solverTypeFromString(std::string &s)
{
	if (s == "NAVIER_STOKES")
		return NAVIER_STOKES;
	else if (s == "SAIKI_BIRINGEN")
		return SAIKI_BIRINGEN;
	else if (s == "FADLUN_ET_AL")
		return FADLUN_ET_AL;
	else if (s == "TAIRA_COLONIUS")
		return TAIRA_COLONIUS;
	else
		return NAVIER_STOKES;
}

SimulationParameters::SimulationParameters(std::string fileName)
{
	PetscInt    rank;
	std::string solver;
	
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	if(rank==0)
	{
		std::ifstream file(fileName.c_str());
		YAML::Parser  parser(file);
		YAML::Node    doc;
		
		parser.GetNextDocument(doc);
				
		doc[0]["dt"] >> dt;
		doc[0]["nt"] >> nt;
		doc[0]["nsave"] >> nsave;
		doc[0]["ibmScheme"] >> solver;
		solverType = solverTypeFromString(solver);
	}
	MPI_Barrier(PETSC_COMM_WORLD);
	
	// broadcast parameters to all processes
	MPI_Bcast(&dt, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
	MPI_Bcast(&nt, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
	MPI_Bcast(&nsave, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
}
