/**
 * \file createKSP.cpp
 * \brief definition of createKSP.
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \date 2017-05-31
 */


// header
# include "createKSP.hpp"


// definition of create KSP
PetscErrorCode createKSP(KSP &ksp, Mat &A, char *FN)
{
    PetscFunctionBeginUser;

    PetscErrorCode  ierr;

    // read KSP settings from a file
    ierr = PetscOptionsInsertFile(PETSC_COMM_WORLD,
            nullptr, FN, PETSC_FALSE); CHKERRQ(ierr);

    // create ksp solver
    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp, A, A); CHKERRQ(ierr);
    ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE); CHKERRQ(ierr);
    ierr = KSPSetType(ksp, KSPCG); CHKERRQ(ierr);
    ierr = KSPSetReusePreconditioner(ksp, PETSC_TRUE); CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
