/*! Implementation of the class `AMGXSolver`.
 * \file amgxsolver.h
 */

#if !defined(AMGXSOLVER_H)
#define AMGXSOLVER_H

#include "solver.h"
#include "AmgXSolver.hpp"

#include <string>


/*!
 * \class AMGXSolver
 * \brief Iterative solver using wrapper for AmgX.
 */
class AMGXSolver : public Solver
{
public:
  AMGXSolver(std::string f): options(f) { };
  virtual ~AMGXSolver(){
    amgx.finalize();
  };

  PetscErrorCode create(const Mat &A);
  PetscErrorCode solve(Vec &x, Vec &b);
  PetscErrorCode getIters(PetscInt &iters);

private:
  AmgXSolver amgx;
  std::string options;

}; // AMGXSolver

#endif