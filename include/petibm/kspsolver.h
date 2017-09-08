/*! Implementation of the class `KSPSolver`.
 * \file kspsolver.h
 */


# pragma once

// PETSc
#include <petscksp.h>

// PetIBM
#include "linsolver.h"


namespace petibm
{
namespace linsolvers
{

/*!
 * \class KSPSolver
 * \brief Iterative solver using PETSc KSP.
 */
class KSPSolver : public LinSolver
{
public:

    /** \copydoc LinSolver::LinSolver(const std::string &, const std::string &). */
    KSPSolver(const std::string &solverName, const std::string &file);

    /** \brief destructor. */
    virtual ~KSPSolver();

    /** \copydoc LinSolver::setMatrix. */
    virtual PetscErrorCode setMatrix(const Mat &A);

    /** \copydoc LinSolver::solve. */
    virtual PetscErrorCode solve(Vec &x, Vec &b);

    /** \copydoc LinSolver::getIters. */
    virtual PetscErrorCode getIters(PetscInt &iters);

    /** \copydoc LinSolver::printInfo. */
    virtual PetscErrorCode printInfo(
            PetscViewer viewer=PETSC_VIEWER_STDOUT_WORLD);

private:

    /** \brief the underlying KSP solver. */
    KSP ksp;

    /** \copydoc LinSolver::init. */
    virtual PetscErrorCode init();

}; // KSPSolver

} // end of namespace linsolvers
} // end of namespace petibm
