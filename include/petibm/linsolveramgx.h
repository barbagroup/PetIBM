/** 
 * \file linsolveramgx.h
 * \brief Def. of LinSolverAmgX.
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \copyright MIT.
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
    
/**
 * \class LinSolverAmgX
 * \brief Iterative solver using wrapper of AmgX.
 * 
 * This class holds a AmgXSolver from AmgXWrapper package. Users must enable
 * AmgX and AmgXWrapper when building PetIBM in order to use this class.
 * 
 * \see petibm::type::LinSolver, petibm::linsolver::createLinSolver.
 * \ingroup linsolver
 */
class LinSolverAmgX : public LinSolverBase
{
public:

    /** \copydoc LinSolverBase(const std::string &, const std::string &) */
    LinSolverAmgX(const std::string &name, const std::string &file);

    /** \copydoc ~LinSolverBase */
    virtual ~LinSolverAmgX();

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

    /** \brief the underlying AmgX wrapper solver. */
    AmgXSolver amgx;

    /** \copydoc LinSolverBase::init */
    virtual PetscErrorCode init();

}; // LinSolverAmgX

} // end of namespace linsolvers
} // end of namespace petibm
