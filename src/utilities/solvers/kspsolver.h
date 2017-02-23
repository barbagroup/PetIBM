/*! Implementation of the class `KSPSolver`.
 * \file kspsolver.h
 */

#if !defined(KSPSOLVER_H)
#define KSPSOLVER_H

#include "solver.h"

#include <petscksp.h>

#include <string>


/*!
 * \class KSPSolver
 * \brief Iterative solver using PETSc KSP.
 */
class KSPSolver : public Solver
{
public:
  KSPSolver(std::string p, std::string f): prefix(p), options(f) { };
  virtual ~KSPSolver(){
    KSPDestroy(&ksp);
  };

  PetscErrorCode create(const Mat &A, const DM &daPack, 
          const PetscBool &split=PETSC_FALSE);
  PetscErrorCode solve(Vec &x, Vec &b);
  PetscErrorCode getIters(PetscInt &iters);

private:
  KSP ksp;
  std::string prefix;
  std::string options;

  PetscErrorCode setSplitPC(const PetscInt &numDM, const DM &daPack);


}; // KSPSolver

#endif
