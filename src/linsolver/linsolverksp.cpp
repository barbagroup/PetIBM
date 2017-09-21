/*! Implementation of the methods of the class `LinSolverKSP`.
 * \file kspsolver.cpp
 */


// PetIBM
#include <petibm/linsolverksp.h>


namespace petibm
{
namespace linsolver
{

// constructor
LinSolverKSP::LinSolverKSP(const std::string &_name, const std::string &_config):
    LinSolverBase(_name, _config) { init(); }


// destructor
LinSolverKSP::~LinSolverKSP() { KSPDestroy(&ksp); }


// underlying initialization function
PetscErrorCode LinSolverKSP::init()
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    ierr = PetscOptionsInsertFile(PETSC_COMM_WORLD, nullptr, 
            config.c_str(), PETSC_FALSE); CHKERRQ(ierr);

    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERRQ(ierr);
    ierr = KSPSetOptionsPrefix(ksp, (name + "_").c_str()); CHKERRQ(ierr);
    ierr = KSPSetType(ksp, KSPCG); CHKERRQ(ierr);
    ierr = KSPSetReusePreconditioner(ksp, PETSC_TRUE); CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


// set coefficient matrix
PetscErrorCode LinSolverKSP::setMatrix(const Mat &A)
{
    PetscFunctionBeginUser;

    PetscErrorCode ierr;

    ierr = KSPSetOperators(ksp, A, A); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // setMatrix


// solve linear system
PetscErrorCode LinSolverKSP::solve(Vec &x, Vec &b)
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


// get number of iterations
PetscErrorCode LinSolverKSP::getIters(PetscInt &iters)
{
    PetscFunctionBeginUser;

    PetscErrorCode ierr;

    ierr = KSPGetIterationNumber(ksp, &iters); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // getIters

} // end of namespace linsolver
} // end of namespace petibm
