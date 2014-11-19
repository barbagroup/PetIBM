/***************************************************************************//**
* \file
* \brief Header to define class SimulationParameters
*/

#if !defined(SIMULATION_PARAMETERS_H)
#define SIMULATION_PARAMETERS_H

#include "types.h"
#include <petscsys.h>
#include <string>

/***************************************************************************//**
* \brief Store various parameters used in the simulation
*/
class SimulationParameters
{
public:
	PetscReal          dt; ///< Size of time increment
	
	PetscInt           nt,        ///< Number of time steps
	                   nsave,     ///< Intervals at which simulation data is saved
	                   startStep; ///< The starting time step of the simulation
	
	SolverType         solverType; ///< Type of flow solver used
	
	TimeSteppingScheme convectionScheme, ///< Time-stepping scheme for the convection term
	                   diffusionScheme;  ///< Time-stepping scheme for the diffusion term
	
	PetscReal          gamma,         ///< Time-stepping coefficient for the convection term in the current time step
	                   zeta,          ///< Time-stepping coefficient for the convection term in the previous time step
	                   alphaExplicit, ///< Time-stepping coefficient for the explicit part of the diffusion term
	                   alphaImplicit; ///< Time-stepping coefficient for the implicit part of the diffusion term
	
	PetscBool          restart; ///< Flag to indicate whether the simulation was restarted from saved data

	PetscReal          velocitySolveTolerance, PoissonSolveTolerance;
	PetscInt           velocitySolveMaxIts, PoissonSolveMaxIts;
	
	/***********************************************************************//**
	* \brief Reads an input file and initializes the simulation parameters
	*/
	SimulationParameters(std::string fileName);
	SimulationParameters();
	void initialize(std::string fileName);
};

#endif