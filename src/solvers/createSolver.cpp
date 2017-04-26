/***************************************************************************//**
 * \file createSolver.cpp
 * \author Anush Krishan (anush@bu.edu)
 * \brief Implementation of the function to create the appropriate solver.
 */


#include "createSolver.h"

/**
 * \brief Creates the appropriate solver.
 *
 * \param cartesianMesh Contains info about the Cartesian grid
 * \param flowDescription Contains info about the fluid and the flow
 * \param simulationParameters Contains info about the parameters of the simulation
 *
 * If there is no immersed boundary in the domain, a Navier-Stokes solver is
 * created. Otherwise, either Taira and Colonius (2007) solver type or
 * Li et al. (2016) solver type is instantiated.
 */
template <PetscInt dim>
std::unique_ptr< NavierStokesSolver<dim> > createSolver(CartesianMesh *cartesianMesh, 
                                                        FlowDescription<dim> *flowDescription, 
                                                        SimulationParameters *simulationParameters)
{
  switch(simulationParameters->ibm)
  {
    case NAVIER_STOKES:
      return std::unique_ptr< NavierStokesSolver<dim> >(new NavierStokesSolver<dim>(cartesianMesh, 
                                                                                    flowDescription, 
                                                                                    simulationParameters));
      break;
    case TAIRA_COLONIUS:
      return std::unique_ptr< TairaColoniusSolver<dim> >(new TairaColoniusSolver<dim>(cartesianMesh, 
                                                                                      flowDescription, 
                                                                                      simulationParameters));
      break;
    case LI_ET_AL:
      return std::unique_ptr< LiEtAlSolver<dim> >(new LiEtAlSolver<dim>(cartesianMesh, 
                                                                        flowDescription, 
                                                                        simulationParameters));
      break;
    default:
      PetscPrintf(PETSC_COMM_WORLD, "ERROR: Unrecognized solver.\n");
      return NULL;
  }
} // createSolver


// dimensions specialization
template std::unique_ptr< NavierStokesSolver<2> > createSolver(CartesianMesh *cartesianMesh, 
                                                               FlowDescription<2> *flowDescription, 
                                                               SimulationParameters *simulationParameters); 
template std::unique_ptr< NavierStokesSolver<3> > createSolver(CartesianMesh *cartesianMesh, 
                                                               FlowDescription<3> *flowDescription, 
                                                               SimulationParameters *simulationParameters);