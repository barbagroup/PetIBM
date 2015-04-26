/***************************************************************************//**
 * \file SimulationParameters.h
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Definition of the class \c SimulationParameters.
 */


#if !defined(SIMULATION_PARAMETERS_H)
#define SIMULATION_PARAMETERS_H

#include "types.h"

#include <string>

#include <petscsys.h>


/**
 * \class SimulationParameters
 * \brief Stores various parameters used in the simulation
 */
class SimulationParameters
{
public:
  PetscReal dt; ///< time-increment
  
  PetscInt startStep, ///< initial time-step 
           nt,        ///< number of time steps
           nsave;     ///< data-saving interval
  
  PetscBool restartFromSolution;  ///< flag to restart from given solution

  SolverType solverType;  ///< type of flow solver
  
  TimeSteppingScheme convectionScheme, ///< time-scheme for the convection term
                     diffusionScheme;  ///< time-scheme for the diffusion term
  
  PetscReal gamma,         ///< coefficient of the convection term at current time step
            zeta,          ///< coefficient of the convection term at previous time step
            alphaExplicit, ///< coefficient of the explicit diffusion term
            alphaImplicit; ///< coefficient of the implicit diffusion term

  PetscReal velocitySolveTolerance, ///< tolerance (velocity solver)
            PoissonSolveTolerance;  ///< tolerance (Poisson solver)
  PetscInt velocitySolveMaxIts, ///< maximum number of iterations (velocity solver)
           PoissonSolveMaxIts;  ///< maximum number of iterations (Poisson solver)
  
  // Parse file and store simulation parameters
  SimulationParameters(std::string fileName);
  SimulationParameters();
  void initialize(std::string fileName);
};

#endif