/***************************************************************************//**
 * \file createSolver.cpp
 * \author Anush Krishan (anush@bu.edu)
 * \brief Implementation of the function to create the appropriate solver.
 */


#include "createSolver.h"

/**
 * \brief Creates the appropriate solver.
 *
 * If there is no immersed boundary in the domain, a Navier-Stokes solver iis
 * created. Otherwise, Taira and Colonius (2007) solver type is instanciated.
 */
template <PetscInt dim>
std::unique_ptr< NavierStokesSolver<dim> > createSolver(std::string directory, 
                                                        CartesianMesh *cartesianMesh, 
                                                        FlowDescription *flowDescription, 
                                                        SimulationParameters *simulationParameters)
{
  if (flowDescription->dimensions != dim)
  {
    PetscPrintf(PETSC_COMM_WORLD, "ERROR: Check the number of dimensions in the input files.\n");
    exit(0);
  }
  switch(simulationParameters->solverType)
  {
    case NAVIER_STOKES:
      return std::unique_ptr< NavierStokesSolver<dim> >(new NavierStokesSolver<dim>(directory,
                                                                                    cartesianMesh, 
                                                                                    flowDescription, 
                                                                                    simulationParameters));
      break;
    case TAIRA_COLONIUS:
      return std::unique_ptr< TairaColoniusSolver<dim> >(new TairaColoniusSolver<dim>(directory,
                                                                                      cartesianMesh, 
                                                                                      flowDescription, 
                                                                                      simulationParameters));
      break;
    default:
      PetscPrintf(PETSC_COMM_WORLD, "ERROR: Unrecognized solver.\n");
      return NULL;
  }
} // createSolver

// specialization
template std::unique_ptr< NavierStokesSolver<2> > createSolver(std::string directory, 
                                                               CartesianMesh *cartesianMesh, 
                                                               FlowDescription *flowDescription, 
                                                               SimulationParameters *simulationParameters); 
template std::unique_ptr< NavierStokesSolver<3> > createSolver(std::string directory, 
                                                               CartesianMesh *cartesianMesh, 
                                                               FlowDescription *flowDescription, 
                                                               SimulationParameters *simulationParameters);