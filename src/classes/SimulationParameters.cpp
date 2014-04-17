#include "SimulationParameters.h"
#include "yaml-cpp/yaml.h"
#include <string>
#include <fstream>

SimulationParameters::SimulationParameters(std::string fileName)
{
	PetscInt       rank;
	
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	// first pass
	if(rank==0)
	{
		std::ifstream file(fileName.c_str());
		YAML::Parser  parser(file);
		YAML::Node    doc;
		
		parser.GetNextDocument(doc);
				
		// cycle through each direction
		for (size_t i=0; i<doc.size(); i++)
		{
			doc[i]["dt"] >> dt;
			doc[i]["nt"] >> nt;
			doc[i]["nsave"] >> nsave;
		}
	}
	MPI_Barrier(PETSC_COMM_WORLD);
	
	// broadcast parameters to all processes
	MPI_Bcast(&dt, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
	MPI_Bcast(&nt, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
	MPI_Bcast(&nsave, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
	
	MPI_Barrier(PETSC_COMM_WORLD);
	
	std::cout << "dt: " << dt << ", nt: " << nt << ", nsave: " << nsave << std::endl;
}
