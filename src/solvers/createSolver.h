/***************************************************************************//**
 * \file createSolver.h
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Definition of the function to create the appropriate solver.
 */


#if !defined(CREATE_SOLVER_H)
#define CREATE_SOLVER_H

#include "navierStokes/NavierStokesSolver.h"
#include "tairaColonius/TairaColoniusSolver.h"

#include <memory>


// create appropriate solver depending on method chosen
template <PetscInt dim>
std::unique_ptr< NavierStokesSolver<dim> > createSolver(CartesianMesh *cartesianMesh, 
                                                        FlowDescription<dim> *flowDescription, 
                                                        SimulationParameters *simulationParameters);

#endif