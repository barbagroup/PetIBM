/**
 * \file io.cpp
 * \brief functions regarding to reading Mat or Vec.
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \date 2017-05-31
 */


// headers
# include "io.hpp"


// definition of readMat
PetscErrorCode readMat(Mat &mat, char *FN, const char *name)
{
    PetscFunctionBeginUser;

    PetscErrorCode  ierr;
    PetscViewer     reader;


    ierr = MatCreate(PETSC_COMM_WORLD, &mat); CHKERRQ(ierr);
    ierr = MatSetType(mat, MATAIJ); CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) mat, name); CHKERRQ(ierr);

    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, 
            FN, FILE_MODE_READ, &reader); CHKERRQ(ierr);

    ierr = MatLoad(mat, reader); CHKERRQ(ierr);

    ierr = PetscViewerDestroy(&reader); CHKERRQ(ierr);

    // get and print matrix information
    {
        MatInfo     info;
        PetscInt    nx,
                    ny;

        ierr = MatGetInfo(mat, MAT_GLOBAL_SUM, &info); CHKERRQ(ierr);
        ierr = MatGetSize(mat, &nx, &ny); CHKERRQ(ierr);

        ierr = PetscPrintf(PETSC_COMM_WORLD, "Matrix %s: %d x %d, with %ld\n", 
                name, nx, ny, (long) info.nz_used); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}


// definition of readVec
PetscErrorCode readVec(Vec &vec, char *FN, const char *name)
{
    PetscFunctionBeginUser;

    PetscErrorCode  ierr;
    PetscViewer     reader;

    ierr = VecCreate(PETSC_COMM_WORLD, &vec); CHKERRQ(ierr);
    ierr = VecSetType(vec, VECSTANDARD); CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) vec, name); CHKERRQ(ierr);

    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, 
            FN, FILE_MODE_READ, &reader); CHKERRQ(ierr);

    ierr = VecLoad(vec, reader); CHKERRQ(ierr);

    ierr = PetscViewerDestroy(&reader); CHKERRQ(ierr);

    // get and print rhs information
    {
        PetscInt    n;
        ierr = VecGetSize(vec, &n); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,
                "Vector %s: size %d\n", name, n); CHKERRQ(ierr);
    }


    PetscFunctionReturn(0);
}


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
