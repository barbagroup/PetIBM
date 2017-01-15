/*! Implementation of the class `Solver`.
 * \file solver.h
 */

#if !defined(SOLVER_H)
#define SOLVER_H

#include <petsc.h>


/**
 * \class Solver
 * \brief Super class for an iterative solver.
 */
class Solver
{
public:
  virtual ~Solver(){ }
  virtual PetscErrorCode create(const Mat &A) = 0;
  virtual PetscErrorCode solve(Vec &x, Vec &b) = 0;
  virtual PetscErrorCode getIters(PetscInt &iters) = 0;

}; // Solver

#endif
