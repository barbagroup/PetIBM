/*! Implementation of the method `createLagrangianSolver` of the class `LiEtAlSolver`.
 * \file createForceSolver.inl
 */

#include "solvers/kspsolver.h"


/*!
 * \brief Creates the solver for the Lagrangian forces.
 */
template <PetscInt dim>
PetscErrorCode LiEtAlSolver<dim>::createForceSolver()
{
  PetscErrorCode ierr;

  std::string prefix;
  std::string options;

  PetscFunctionBeginUser;

  // possibility to overwrite the path of the configuration file
  // using the command-line parameter: `-forces_config_file <filepath>`
  char path[PETSC_MAX_PATH_LEN];
  PetscBool found;
  PetscOptionsGetString(NULL, NULL, "-forces_config_file", path, sizeof(path), &found);

  prefix = "forces_";
  options = (found) ? std::string(path) : NavierStokesSolver<dim>::parameters->directory + "/solversPetscOptions.info";
  forces = new KSPSolver(prefix, options);
  forces->create(EBNET);

  PetscFunctionReturn(0);
} // createForceSolver
