/***************************************************************************//**
 * \file DiffusionTerm.h
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \brief Definition of the class \c DiffusionTerm.
 */


 #if !defined(DIFFUSION_TERM_H)
 #define DIFFUSION_TERM_H

 #include <NavierStokesSolver.h>


 /**
  * \class DiffusionTerm
  * \brief Computes the numerical and exact explicit diffusion terms.
  */
template <PetscInt dim>
class DiffusionTerm : public NavierStokesSolver<dim>
{
public:
  Vec rnExact;  // exact solution of the explicit diffusion term
  PetscReal relativeError;  // relative error in the explicit diffusion term

  PetscErrorCode initialize();
  PetscErrorCode initializeFluxes();
  PetscErrorCode calculateExactSolution();
  PetscErrorCode calculateRelativeError();
  PetscErrorCode writeRelativeError();
  PetscErrorCode finalize();

  DiffusionTerm(std::string folder, FlowDescription *FD, SimulationParameters *SP, CartesianMesh *CM);

  virtual std::string name()
  {
    return "Difffusion term";
  }
};

#endif