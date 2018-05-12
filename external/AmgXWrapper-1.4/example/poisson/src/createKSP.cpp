/**
 * \file createKSP.cpp
 * \brief prototype of createKSP
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \date 2017-05-31
 */


// header
# include "createKSP.hpp"


// definition of createKSP
PetscErrorCode createKSP(KSP &ksp, Mat &A, DM &grid, char *FN)
{
    PetscFunctionBeginUser;

    PetscErrorCode  ierr;

    // create ksp solver
    ierr = PetscOptionsInsertFile(PETSC_COMM_WORLD,
            nullptr, FN, PETSC_FALSE); CHKERRQ(ierr);

    // create KSP for intermaediate fluxes system
    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERRQ(ierr);
    ierr = KSPSetDM(ksp, grid); CHKERRQ(ierr);
    ierr = KSPSetDMActive(ksp, PETSC_FALSE); CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp, A, A); CHKERRQ(ierr);
    ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE); CHKERRQ(ierr);
    ierr = KSPSetType(ksp, KSPCG); CHKERRQ(ierr);
    ierr = KSPSetReusePreconditioner(ksp, PETSC_TRUE); CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
