/*! Implementation of the methods of the class `AMGXSolver`.
 * \file amgxsolver.cpp
 */

#include "amgxsolver.h"


AMGXSolver::~AMGXSolver()
{
	PetscErrorCode ierr;
	PetscBool finalized;

	PetscFunctionBeginUser;

	ierr = PetscFinalized(&finalized); CHKERRV(ierr);
	if (finalized) return;

	ierr = amgx.finalize(); CHKERRV(ierr);
} // ~AMGXSolver


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
  PetscErrorCode ierr;

  ierr = amgx.getIters(iters); CHKERRQ(ierr);

  return 0;
} // getIters
