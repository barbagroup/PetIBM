/***************************************************************************//**
 * \file createSolver.h
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Definition of the function to create the appropriate solver.
 */


#if !defined(CREATE_SOLVER_H)
#define CREATE_SOLVER_H

#include "NavierStokesSolver.h"
#include "TairaColoniusSolver.h"

#include <memory>


template <PetscInt dim>
std::unique_ptr< NavierStokesSolver<dim> > createSolver(std::string folder, 
                                                        FlowDescription *FD, 
                                                        SimulationParameters *SP, 
                                                        CartesianMesh *CM);

#endif