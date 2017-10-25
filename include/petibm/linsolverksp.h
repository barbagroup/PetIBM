/*! Implementation of the class `KSPSolver`.
 * \file kspsolver.h
 */


# pragma once

// PETSc
# include <petscksp.h>

// PetIBM
# include <petibm/linsolver.h>


namespace petibm
{
namespace linsolver
{
/*!
 * \class LinSolverKSP
 * \brief Iterative solver using PETSc KSP.
 */
class LinSolverKSP : public LinSolverBase
{
public:

    /** \copydoc LinSolverBase::LinSolverBase. */
    LinSolverKSP(const std::string &name, const std::string &file);

    /** \brief destructor. */
    virtual ~LinSolverKSP();

    /** \copydoc LinSolverBase::setMatrix. */
    virtual PetscErrorCode setMatrix(const Mat &A);

    /** \copydoc LinSolverBase::solve. */
    virtual PetscErrorCode solve(Vec &x, Vec &b);

    /** \copydoc LinSolverBase::getIters. */
    virtual PetscErrorCode getIters(PetscInt &iters);

    /** \copydoc LinSolverBase::getResidual. */
    virtual PetscErrorCode getResidual(PetscReal &res);

private:

    /** \brief the underlying KSP solver. */
    KSP ksp;

    /** \copydoc LinSolverBase::init. */
    virtual PetscErrorCode init();

}; // LinSolverKSP

} // end of namespace linsolver
} // end of namespace petibm
