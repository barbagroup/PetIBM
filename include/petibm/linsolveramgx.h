/*! Implementation of the class `LinSolverAmgX`.
 * \file amgxsolver.h
 */


# pragma once

// AmgX wrapper
# include <AmgXSolver.hpp>

// PetIBM
# include <petibm/linsolver.h>


namespace petibm
{
namespace linsolver
{
/*!
 * \class LinSolverAmgX
 * \brief Iterative solver using wrapper for AmgX.
 */
class LinSolverAmgX : public LinSolverBase
{
public:

    /** \copydoc LinSolverBase::LinSolverBase. */
    LinSolverAmgX(const std::string &name, const std::string &file);

    /** \brief destructor. */
    virtual ~LinSolverAmgX();

    /** \copydoc LinSolverBase::setMatrix. */
    virtual PetscErrorCode setMatrix(const Mat &A);

    /** \copydoc LinSolverBase::solve. */
    virtual PetscErrorCode solve(Vec &x, Vec &b);

    /** \copydoc LinSolverBase::getIters. */
    virtual PetscErrorCode getIters(PetscInt &iters);

    /** \copydoc LinSolverBase::getResidual. */
    virtual PetscErrorCode getResidual(PetscReal &res);

private:

    /** \brief the underlying AmgX wrapper. */
    AmgXSolver amgx;

    /** \copydoc LinSolverBase::init. */
    virtual PetscErrorCode init();

}; // LinSolverAmgX

} // end of namespace linsolvers
} // end of namespace petibm
