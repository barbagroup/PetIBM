/**
 * \file io.cpp
 * \brief functions regarding to reading Mat or Vec.
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \date 2017-05-31
 */


// headers
# include "io.hpp"


// definition of printHeader
PetscErrorCode printHeader(StructArgs &args)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    ierr = PetscPrintf(PETSC_COMM_WORLD, "\n"); CHKERRQ(ierr);
    for(int i=0; i<72; ++i) ierr = PetscPrintf(PETSC_COMM_WORLD, "=");
    ierr = PetscPrintf(PETSC_COMM_WORLD, "\n"); CHKERRQ(ierr);
    for(int i=0; i<72; ++i) ierr = PetscPrintf(PETSC_COMM_WORLD, "-");
    ierr = PetscPrintf(PETSC_COMM_WORLD, "\n"); CHKERRQ(ierr);

    ierr = args.print(); CHKERRQ(ierr);

    for(int i=0; i<72; ++i) ierr = PetscPrintf(PETSC_COMM_WORLD, "-");
    ierr = PetscPrintf(PETSC_COMM_WORLD, "\n"); CHKERRQ(ierr);
    for(int i=0; i<72; ++i) ierr = PetscPrintf(PETSC_COMM_WORLD, "=");
    ierr = PetscPrintf(PETSC_COMM_WORLD, "\n"); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


// definition of printHeader
PetscErrorCode printFooter(StructArgs &args)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    ierr = PetscPrintf(PETSC_COMM_WORLD, "\n"); CHKERRQ(ierr);
    for(int i=0; i<72; ++i) ierr = PetscPrintf(PETSC_COMM_WORLD, "=");
    ierr = PetscPrintf(PETSC_COMM_WORLD, "\n"); CHKERRQ(ierr);
    for(int i=0; i<72; ++i) ierr = PetscPrintf(PETSC_COMM_WORLD, "-");
    ierr = PetscPrintf(PETSC_COMM_WORLD, "\n"); CHKERRQ(ierr);

    ierr = PetscPrintf(PETSC_COMM_WORLD, "End of %s\n", args.caseName); CHKERRQ(ierr);

    for(int i=0; i<72; ++i) ierr = PetscPrintf(PETSC_COMM_WORLD, "-");
    ierr = PetscPrintf(PETSC_COMM_WORLD, "\n"); CHKERRQ(ierr);
    for(int i=0; i<72; ++i) ierr = PetscPrintf(PETSC_COMM_WORLD, "=");
    ierr = PetscPrintf(PETSC_COMM_WORLD, "\n"); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
