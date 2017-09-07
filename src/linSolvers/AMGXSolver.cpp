/*! Implementation of the methods of the class `AMGXSolver`.
 * \file amgxsolver.cpp
 */


// PetIBM
#include "AMGXSolver.h"


namespace petibm
{
namespace linsolvers
{

/** \copydoc KSPSolver::KSPSolver. */
AMGXSolver::AMGXSolver(const std::string &_name, const std::string &_options):
    LinSolver(_name, _options) { init(); }


/** \copydoc KSPSolver::~KSPSolver. */
AMGXSolver::~AMGXSolver() { amgx.finalize(); }


/** \copydoc AMGXSolver::init. */
PetscErrorCode AMGXSolver::init()
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    ierr = amgx.initialize(PETSC_COMM_WORLD, "dDDI", options); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/** \copydoc AMGXSolver::setMatrix. */
PetscErrorCode AMGXSolver::setMatrix(const Mat &A)
{
    PetscErrorCode ierr;

    ierr = amgx.setA(A); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // setMatrix


/** \copydoc AMGXSolver::solve. */
PetscErrorCode AMGXSolver::solve(Vec &x, Vec &b)
{
    PetscErrorCode ierr;

    ierr = amgx.solve(x, b); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // solve


/** \copydoc AMGXSolver::getIters. */
PetscErrorCode AMGXSolver::getIters(PetscInt &iters)
{
    PetscErrorCode ierr;

    iters = amgx.getIters();

    PetscFunctionReturn(0);
} // getIters

/** \copydoc AMGXSolver::printInfo. */
PetscErrorCode AMGXSolver::printInfo(PetscViewer viewer)
{
    PetscFunctionBeginUser;
    // TODO: implement it.
    PetscFunctionReturn(0);
}

} // end of namespace linsolvers
} // end of namespace petibm
