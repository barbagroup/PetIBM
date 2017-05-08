/*! Implementation of the class `AMGXSolver`.
 * \file amgxsolver.h
 */


# pragma once
#define AMGXSOLVER_H


// AmgX wrapper
# include "AmgXSolver.hpp"

// PetIBM
# include "LinSolver.h"


/*!
 * \class AMGXSolver
 * \brief Iterative solver using wrapper for AmgX.
 */
class AMGXSolver : public LinSolver
{
public:

    /** \copydoc LinSolver::LinSolver(const std::string &, const std::string &). */
    AMGXSolver(const std::string &solverName, const std::string &file);

    /** \brief destructor. */
    virtual ~AMGXSolver();

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

    /** \brief the underlying AmgX wrapper. */
    AmgXSolver amgx;

    /** \copydoc LinSolver::init. */
    virtual PetscErrorCode init();

}; // AMGXSolver
