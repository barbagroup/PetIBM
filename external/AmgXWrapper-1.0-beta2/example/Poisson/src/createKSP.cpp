# include "headers.hpp"


PetscErrorCode createKSP(KSP &ksp, Mat &A, DM &grid, char *FN)
{
    PetscErrorCode  ierr;

    // create ksp solver
    ierr = PetscOptionsInsertFile(
            PETSC_COMM_WORLD, nullptr, FN, PETSC_FALSE);                CHKERRQ(ierr);

    // create KSP for intermaediate fluxes system
    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);                           CHKERRQ(ierr);
    ierr = KSPSetDM(ksp, grid);                                         CHKERRQ(ierr);
    ierr = KSPSetDMActive(ksp, PETSC_FALSE);                            CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp, A, A);                                  CHKERRQ(ierr);
    ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);                  CHKERRQ(ierr);
    ierr = KSPSetType(ksp, KSPCG);                                      CHKERRQ(ierr);
    ierr = KSPSetReusePreconditioner(ksp, PETSC_TRUE);                  CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp);                                      CHKERRQ(ierr);

    return 0;
}
