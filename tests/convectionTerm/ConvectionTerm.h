/***************************************************************************//**
 * \file convectionTerm.h
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \brief Definition of the class \c ConvectionTerm.
 */


 #if !defined(CONVECTION_TERM_H)
 #define CONVECTION_TERM_H

 #include <NavierStokesSolver.h>


 /**
  * \class ConvectionTerm
  * \brief Computes the numerical and exact explicit diffusion terms.
  */
template <PetscInt dim>
class ConvectionTerm : public NavierStokesSolver<dim>
{
public:
  Vec HExact;               // exact solution of the explicit convection term
  PetscReal relativeError;  // relative error in the explicit convection term

  PetscErrorCode initialize();
  PetscErrorCode initializeFluxes();
  PetscErrorCode calculateExactSolution();
  PetscErrorCode calculateRelativeError();
  PetscErrorCode writeRelativeError();
  PetscErrorCode finalize();

  ConvectionTerm(std::string folder, FlowDescription *FD, SimulationParameters *SP, CartesianMesh *CM);

  virtual std::string name()
  {
    return "Convection term";
  }
};

#endif