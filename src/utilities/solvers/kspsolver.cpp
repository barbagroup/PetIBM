/*! Implementation of the methods of the class `KSPSolver`.
 * \file kspsolver.cpp
 */

#include "kspsolver.h"


/*!
 * \brief Creates the KSP solver.
 */
PetscErrorCode KSPSolver::create(const Mat &A)
{
  PetscErrorCode ierr;

  ierr = PetscOptionsInsertFile(PETSC_COMM_WORLD, NULL, 
                                options.c_str(), PETSC_FALSE); CHKERRQ(ierr);

  // create KSP for intermediate fluxes system
  ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERRQ(ierr);
  ierr = KSPSetOptionsPrefix(ksp, prefix.c_str()); CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp, A, A); CHKERRQ(ierr);
  ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE); CHKERRQ(ierr);
  ierr = KSPSetType(ksp, KSPCG); CHKERRQ(ierr);
  ierr = KSPSetReusePreconditioner(ksp, PETSC_TRUE); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);

  return 0;
} // create


/*!
 * \brief Solves the system.
 */
PetscErrorCode KSPSolver::solve(Vec &x, Vec &b)
{
  PetscErrorCode ierr;

  ierr = KSPSolve(ksp, b, x); CHKERRQ(ierr);

  KSPConvergedReason reason;
  ierr = KSPGetConvergedReason(ksp, &reason); CHKERRQ(ierr);
  if (reason < 0)
  {
    ierr = PetscPrintf(PETSC_COMM_WORLD, 
            "\nERROR: %s solver diverged due to reason: %d\n", 
            prefix.substr(0, prefix.size()-1).c_str(), reason); CHKERRQ(ierr);
    ierr = PetscFinalize(); CHKERRQ(ierr);
    exit(1);
  }

  return 0;
} // solve


/*!
 * \brief Get the number of iterations performed.
 */
PetscErrorCode KSPSolver::getIters(PetscInt &iters)
{
  PetscErrorCode ierr;

  ierr = KSPGetIterationNumber(ksp, &iters); CHKERRQ(ierr);

  return 0;
} // getIters
