/*! Implementation of the method `createSolvers` of the class `NavierStokesSolver.`
 * \file createSolvers.inl
 */

#include "solvers/kspsolver.h"
#ifdef HAVE_AMGX
#include "solvers/amgxsolver.h"
#endif


/*!
 * \brief Creates the iterative solvers.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::createSolvers()
{
  PetscErrorCode ierr;

  std::string prefix;
  std::string options;
  PetscInt idx;

  // possibility to overwrite the path of the configuration file
  // using the command-line parameter: `-velocity_config_file <file-path>`
  char path[PETSC_MAX_PATH_LEN];
  PetscBool found;
  PetscOptionsGetString(NULL, NULL, "-velocity_config_file", path, sizeof(path), &found);

  switch(parameters->vSolveType)
  {
    case CPU:
      prefix = "velocity_";
      options = (found) ? std::string(path) : parameters->directory + "/solversPetscOptions.info";
      velocity = new KSPSolver(prefix, options);
      break;
    case GPU:
#ifdef HAVE_AMGX
      options = (found) ? std::string(path) : parameters->directory + "/solversAmgXOptions_v.info";
      velocity = new AMGXSolver(options);
#else HAVE_AMGX
      ierr = PetscPrintf(PETSC_COMM_WORLD,
                         "\nERROR: AmgX not available; vSolveType should be set to 'CPU'.\n");
      exit(1);
#endif
      break;
    default:
      ierr = PetscPrintf(PETSC_COMM_WORLD,
                         "\nERROR: vSolveType should be either 'CPU' or 'GPU'.\n");
      exit(1);
  }
  velocity->create(A);

  // possibility to overwrite the path of the configuration file
  // using the command-line parameter: `-poisson_config_file <file-path>`
  PetscOptionsGetString(NULL, NULL, "-poisson_config_file", path, sizeof(path), &found);

  switch(parameters->pSolveType)
  {
    case CPU:
      prefix = "poisson_";
      options = (found) ? std::string(path) : parameters->directory + "/solversPetscOptions.info";
      poisson = new KSPSolver(prefix, options);
      ierr = setNullSpace(); CHKERRQ(ierr);
      break;
    case GPU:
#ifdef HAVE_AMGX
      options = (found) ? std::string(path) : parameters->directory + "/solversAmgXOptions_p.info";
      poisson = new AMGXSolver(options);
      // Set the first row of the Poisson matrix to zero (except diagonal element)
      // to fix the pressure at the bottom-left corner of the domain
      // This is done to deal with the null space of the Poisson matrix when using
      // AmgX solver
      idx = 0;
      ierr = MatZeroRows(QTBNQ, 1, &idx, 1.0, PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);
      ierr = PetscObjectViewFromOptions((PetscObject) QTBNQ, NULL, "-QTBNQ_mat_view"); CHKERRQ(ierr);
#else HAVE_AMGX
      ierr = PetscPrintf(PETSC_COMM_WORLD,
                         "\nERROR: AmgX not available; pSolveType should be set to 'CPU'.\n");
      exit(1);
#endif
      break;
    default:
      ierr = PetscPrintf(PETSC_COMM_WORLD,
                         "\nERROR: pSolveType should be either 'CPU' or 'GPU'.\n");
      exit(1);
  }
  poisson->create(QTBNQ);

  return 0;
} // createSolvers
