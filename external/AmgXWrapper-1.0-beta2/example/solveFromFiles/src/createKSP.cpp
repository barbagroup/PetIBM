# include "headers.hpp"


PetscErrorCode createKSP(KSP &ksp, Mat &A, char *FN)
{
    PetscErrorCode  ierr;

    // create ksp solver
    ierr = PetscOptionsInsertFile(
            PETSC_COMM_WORLD, nullptr, FN, PETSC_FALSE);                     CHK;

    // create KSP for intermaediate fluxes system
    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);                                CHK;
    ierr = KSPSetOperators(ksp, A, A);                                       CHK;
    ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);                       CHK;
    ierr = KSPSetType(ksp, KSPCG);                                           CHK;
    ierr = KSPSetReusePreconditioner(ksp, PETSC_TRUE);                       CHK;
    ierr = KSPSetFromOptions(ksp);                                           CHK;

    return 0;
}
