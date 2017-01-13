/*! Implementation of the method `createSolvers` of the class `NavierStokesSolver.`
 * \file createSolvers.inl
 */

#include "solvers/kspsolver.h"
#ifdef HAVE_AMGX
#include "solvers/amgxsolver.h"
#endif


template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::createSolvers()
{
  PetscErrorCode ierr;

  std::string prefix;
  std::string options;

  switch(parameters->vSolveType)
  {
    case CPU:
      prefix = "velocity_";
      options = parameters->directory + "/solversPetscOptions.info";
      velocity = new KSPSolver(prefix, options);
      break;
    case GPU:
#ifdef HAVE_AMGX
      options = parameters->directory + "/solversAmgXOptions_v.info";
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

  switch(parameters->pSolveType)
  {
    case CPU:
      prefix = "poisson_";
      options = parameters->directory + "/solversPetscOptions.info";
      poisson = new KSPSolver(prefix, options);
      break;
    case GPU:
#ifdef HAVE_AMGX
      options = parameters->directory + "/solversAmgXOptions_p.info";
      poisson = new AMGXSolver(options);
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
