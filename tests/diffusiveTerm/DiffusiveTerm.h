/***************************************************************************//**
 * \file DiffusiveTerm.h
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \brief Definition of the class \c DiffusiveTerm.
 */


#if !defined(DIFFUSIVE_TERM_H)
#define DIFFUSIVE_TERM_H

#include <navierStokes/NavierStokesSolver.h>


/**
 * \class DiffusiveTerm
 * \brief Computes the numerical and exact explicit diffusion terms.
 */
template <PetscInt dim>
class DiffusiveTerm : public NavierStokesSolver<dim>
{
public:
  Vec rnExact;  // exact solution of the explicit diffusion term
  PetscReal relativeError;  // relative error in the explicit diffusion term

  PetscErrorCode initializeFluxes();
  PetscErrorCode calculateExactSolution();
  PetscErrorCode calculateRelativeError();
  PetscErrorCode writeRelativeError();
  PetscErrorCode finalize();

  DiffusiveTerm(CartesianMesh *cartesianMesh, 
                FlowDescription<dim> *flowDescription, 
                SimulationParameters *simulationParameters);

}; // DiffusiveTerm

#endif