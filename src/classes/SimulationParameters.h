#if !defined(SIMULATION_PARAMETERS_H)
#define SIMULATION_PARAMETERS_H

#include <petscsys.h>
#include <string>

class SimulationParameters
{
public:
	PetscReal dt;
	PetscInt  nt, nsave;
	
	SimulationParameters(std::string fileName);
};

#endif
