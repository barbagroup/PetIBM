/**
 * \file linsolverksp.h
 * \brief Implementation of LinSolverKSP.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 */

// PetIBM
#include <petibm/linsolverksp.h>

namespace petibm
{
namespace linsolver
{
// implement LinSolverKSP::LinSolverKSP
LinSolverKSP::LinSolverKSP(const std::string &_name, const std::string &_config)
    : LinSolverBase(_name, _config)
{
    init();
}  // LinSolverKSP

// implement LinSolverKSP::~LinSolverKSP
LinSolverKSP::~LinSolverKSP()
{
    PetscFunctionBeginUser;

    PetscErrorCode ierr;
    PetscBool finalized;

    ierr = PetscFinalized(&finalized); CHKERRV(ierr);
    if (finalized) return;

    ierr = KSPDestroy(&ksp); CHKERRV(ierr);
}  // ~LinSolverKSP

// implement LinSolverKSP::destroy
PetscErrorCode LinSolverKSP::destroy()
{
    PetscErrorCode ierr;

    ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
    ierr = LinSolverBase::destroy(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // destroy

// implement LinSolverKSP::init
PetscErrorCode LinSolverKSP::init()
{
    PetscFunctionBeginUser;

    PetscErrorCode ierr;

    type = "PETSc KSP";

    if (config != "None")
    {
        ierr = PetscOptionsInsertFile(PETSC_COMM_WORLD, nullptr, config.c_str(),
                                      PETSC_TRUE); CHKERRQ(ierr);
    }

    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERRQ(ierr);
    ierr = KSPSetOptionsPrefix(ksp, (name + "_").c_str()); CHKERRQ(ierr);
    ierr = KSPSetType(ksp, KSPCG); CHKERRQ(ierr);
    ierr = KSPSetReusePreconditioner(ksp, PETSC_TRUE); CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // init

// implement LinSolverKSP::setMatrix
PetscErrorCode LinSolverKSP::setMatrix(const Mat &A)
{
    PetscFunctionBeginUser;

    PetscErrorCode ierr;

    ierr = KSPReset(ksp); CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp, A, A); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // setMatrix

// implement LinSolverKSP::solve
PetscErrorCode LinSolverKSP::solve(Vec &x, Vec &b)
{
    PetscFunctionBeginUser;

    PetscErrorCode ierr;
    KSPConvergedReason reason;

    ierr = KSPSolve(ksp, b, x); CHKERRQ(ierr);

    ierr = KSPGetConvergedReason(ksp, &reason); CHKERRQ(ierr);

    if (reason < 0)
    {
        ierr = KSPReasonView(ksp, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

        SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_CONV_FAILED,
                 "PetIBM exited due to PETSc KSP solver %s diverged with "
                 "reason %d.",
                 name.c_str(), reason);
    }

    PetscFunctionReturn(0);
}  // solve

// implement LinSolverKSP::getIters
PetscErrorCode LinSolverKSP::getIters(PetscInt &iters)
{
    PetscFunctionBeginUser;

    PetscErrorCode ierr;

    ierr = KSPGetIterationNumber(ksp, &iters); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // getIters

// implement LinSolverKSP::getResidual
PetscErrorCode LinSolverKSP::getResidual(PetscReal &res)
{
    PetscFunctionBeginUser;

    PetscErrorCode ierr;

    ierr = KSPGetResidualNorm(ksp, &res); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // getResidual

}  // end of namespace linsolver
}  // end of namespace petibm
