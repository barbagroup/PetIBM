/*! Implementation of the methods of the class `KSPSolver`.
 * \file kspsolver.cpp
 */


// PetIBM
#include "petibm/kspsolver.h"


namespace petibm
{
namespace linsolvers
{

/** \copydoc KSPSolver::KSPSolver. */
KSPSolver::KSPSolver(const std::string &_name, const std::string &_options):
    LinSolver(_name, _options) { init(); }


/** \copydoc KSPSolver::~KSPSolver. */
KSPSolver::~KSPSolver() { KSPDestroy(&ksp); }


/** \copydoc KSPSolver::init. */
PetscErrorCode KSPSolver::init()
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    ierr = PetscOptionsInsertFile(PETSC_COMM_WORLD, nullptr, 
            options.c_str(), PETSC_FALSE); CHKERRQ(ierr);

    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERRQ(ierr);
    ierr = KSPSetOptionsPrefix(ksp, (name + "_").c_str()); CHKERRQ(ierr);
    ierr = KSPSetType(ksp, KSPCG); CHKERRQ(ierr);
    ierr = KSPSetReusePreconditioner(ksp, PETSC_TRUE); CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/** \copydoc KSPSolver::setMatrix. */
PetscErrorCode KSPSolver::setMatrix(const Mat &A)
{
    PetscFunctionBeginUser;

    PetscErrorCode ierr;

    ierr = KSPSetOperators(ksp, A, A); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // setMatrix


/** \copydoc KSPSolver::solve. */
PetscErrorCode KSPSolver::solve(Vec &x, Vec &b)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;
    KSPConvergedReason  reason;

    ierr = KSPSolve(ksp, b, x); CHKERRQ(ierr);

    ierr = KSPGetConvergedReason(ksp, &reason); CHKERRQ(ierr);

    if (reason < 0)
    {
        ierr = KSPReasonView(ksp, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

        SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_CONV_FAILED,
                "PetIBM exited due to PETSc KSP solver %s dirveged with "
                "reason %d.", name.c_str(), reason);
    }

    return 0;
} // solve


/** \copydoc KSPSolver::getIters. */
PetscErrorCode KSPSolver::getIters(PetscInt &iters)
{
    PetscFunctionBeginUser;

    PetscErrorCode ierr;

    ierr = KSPGetIterationNumber(ksp, &iters); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // getIters


/** \copydoc KSPSolver::printInfo. */
PetscErrorCode KSPSolver::printInfo(PetscViewer viewer)
{
    PetscFunctionBeginUser;

    PetscErrorCode  ierr;

    ierr = KSPView(ksp, viewer); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

} // end of namespace linsolvers
} // end of namespace petibm
