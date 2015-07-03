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
std::unique_ptr< NavierStokesSolver<dim> > createSolver(std::string folder, 
                                                        FlowDescription *FD, 
                                                        SimulationParameters *SP, 
                                                        CartesianMesh *CM)
{
  if(FD->dimensions!=dim)
  {
    PetscPrintf(PETSC_COMM_WORLD, "ERROR: Check the number of dimensions in the input files.\n");
    exit(0);
  }
  switch(SP->solverType)
  {
    case NAVIER_STOKES:
      return std::unique_ptr< NavierStokesSolver<dim> >(new NavierStokesSolver<dim>(folder, 
                                                                                    FD, SP, CM));
      break;
    case TAIRA_COLONIUS:
      return std::unique_ptr< TairaColoniusSolver<dim> >(new TairaColoniusSolver<dim>(folder, 
                                                                                      FD, SP, CM));
      break;
    default:
      PetscPrintf(PETSC_COMM_WORLD, "Unrecognized solver!\n");
      return NULL;
  }
} // createSolver

template std::unique_ptr< NavierStokesSolver<2> > createSolver(std::string folder, 
                                                               FlowDescription *FD, 
                                                               SimulationParameters *SP, 
                                                               CartesianMesh *CM);
template std::unique_ptr< NavierStokesSolver<3> > createSolver(std::string folder, 
                                                               FlowDescription *FD, 
                                                               SimulationParameters *SP, 
                                                               CartesianMesh *CM);