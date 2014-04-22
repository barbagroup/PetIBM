#if !defined(SIMULATION_PARAMETERS_H)
#define SIMULATION_PARAMETERS_H

#include "types.h"
#include <petscsys.h>
#include <string>

class SimulationParameters
{
public:
	PetscReal  dt;
	PetscInt   nt, nsave;
	SolverType solverType;
	SimulationParameters(std::string fileName);
};

#endif
