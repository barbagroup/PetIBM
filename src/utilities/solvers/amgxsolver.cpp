/*! Implementation of the methods of the class `AMGXSolver`.
 * \file amgxsolver.cpp
 */

#include "amgxsolver.h"


/*!
 * \brief Creates the AmgX solver.
 */
PetscErrorCode AMGXSolver::create(const Mat &A)
{
  PetscErrorCode ierr;

  ierr = amgx.initialize(PETSC_COMM_WORLD, "dDDI", options); CHKERRQ(ierr);
  ierr = amgx.setA(A); CHKERRQ(ierr);

  return 0;
} // create


/*!
 * \brief Solves the system.
 */
PetscErrorCode AMGXSolver::solve(Vec &x, Vec &b)
{
  PetscErrorCode ierr;

  ierr = amgx.solve(x, b); CHKERRQ(ierr);

  return 0;
} // solve


/*!
 * \brief Get the number of iterations performed.
 */
PetscErrorCode AMGXSolver::getIters(PetscInt &iters)
{
  iters = amgx.getIters();

  return 0;
} // getIters
