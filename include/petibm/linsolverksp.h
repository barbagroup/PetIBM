/**
 * \file linsolverksp.h
 * \brief Def. of LinSolverKSP.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 */

#pragma once

#include <petscksp.h>

#include <petibm/linsolver.h>

namespace petibm
{
namespace linsolver
{
/**
 * \class LinSolverKSP
 * \brief Iterative solver using PETSc KSP.
 *
 * This class holds a KSP object from PETSc. The configuration of this
 * underlying KSP solver will be read from the provided argument `file`.
 *
 * \see petibm::type::LinSolver, petibm::linsolver::createLinSolver.
 * \ingroup linsolver
 */
class LinSolverKSP : public LinSolverBase
{
public:
    /** \copydoc LinSolverBase(const std::string &, const std::string &)
     *
     * The argument `name` will be used as a prefix for the configuration of the
     * underlying KSP solver in the provided configuration file.
     */
    LinSolverKSP(const std::string &solverName, const std::string &file);

    /** \copydoc ~LinSolverBase */
    virtual ~LinSolverKSP();

    /** \copydoc LinSolverBase::destroy */
    virtual PetscErrorCode destroy();

    /** \copydoc LinSolverBase::setMatrix */
    virtual PetscErrorCode setMatrix(const Mat &A);

    /** \copydoc LinSolverBase::solve */
    virtual PetscErrorCode solve(Vec &x, Vec &b);

    /** \copydoc LinSolverBase::getIters */
    virtual PetscErrorCode getIters(PetscInt &iters);

    /** \copydoc LinSolverBase::getResidual */
    virtual PetscErrorCode getResidual(PetscReal &res);

protected:
    /** \brief the underlying KSP solver */
    KSP ksp;

    /** \copydoc LinSolverBase::init */
    virtual PetscErrorCode init();

};  // LinSolverKSP

}  // end of namespace linsolver

}  // end of namespace petibm
