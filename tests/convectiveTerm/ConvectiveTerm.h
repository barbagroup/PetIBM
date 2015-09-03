/***************************************************************************//**
 * \file ConvectiveTerm.h
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \brief Definition of the class \c ConvectiveTerm.
 */


#if !defined(CONVECTIVE_TERM_H)
#define CONVECTIVE_TERM_H

#include <navierStokes/NavierStokesSolver.h>


/**
 * \class ConvectiveTerm
 * \brief Computes the numerical and exact explicit convective terms.
 */
template <PetscInt dim>
class ConvectiveTerm : public NavierStokesSolver<dim>
{
public:
  Vec rnExact;               // exact solution of the explicit convective term
  PetscReal relativeError;   // relative error in the explicit convective term

  PetscErrorCode initializeFluxes();
  PetscErrorCode calculateExactSolution();
  PetscErrorCode calculateRelativeError();
  PetscErrorCode writeRelativeError();
  PetscErrorCode finalize();

  ConvectiveTerm(CartesianMesh *cartesianMesh,  
                 FlowDescription<dim> *flowDescription, 
                 SimulationParameters *simulationParameters);

}; // ConvectiveTerm

#endif