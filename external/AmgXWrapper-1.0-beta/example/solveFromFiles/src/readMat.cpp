# include "headers.hpp"

PetscErrorCode readMat(Mat &mat, char *FN, const char *name)
{
    PetscErrorCode  ierr;
    PetscViewer     reader;


    ierr = MatCreate(PETSC_COMM_WORLD, &mat);                                CHK;
    ierr = MatSetType(mat, MATAIJ);                                          CHK;
    //ierr = MatSetOption(mat, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);      CHK;
    //ierr = MatSetOption(mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);    CHK;
    ierr = PetscObjectSetName((PetscObject) mat, name);                      CHK;

    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, 
            FN, FILE_MODE_READ, &reader);                                    CHK;

    ierr = MatLoad(mat, reader);                                             CHK;

    ierr = PetscViewerDestroy(&reader);                                      CHK;

    // get and print matrix information
    {
        MatInfo     info;
        PetscInt    nx,
                    ny;

        ierr = MatGetInfo(mat, MAT_GLOBAL_SUM, &info);                       CHK;
        ierr = MatGetSize(mat, &nx, &ny);                                    CHK;

        ierr = PetscPrintf(PETSC_COMM_WORLD, 
                "Matrix %s: %d x %d, with %ld\n", 
                name, nx, ny, (long) info.nz_used);                          CHK;
    }

    return 0;
}
