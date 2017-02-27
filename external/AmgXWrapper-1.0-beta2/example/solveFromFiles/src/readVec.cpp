# include "headers.hpp"

PetscErrorCode readVec(Vec &vec, char *FN, const char *name)
{
    PetscErrorCode  ierr;
    PetscViewer     reader;

    ierr = VecCreate(PETSC_COMM_WORLD, &vec);                                CHK;
    ierr = VecSetType(vec, VECSTANDARD);                                     CHK;
    ierr = PetscObjectSetName((PetscObject) vec, name);                      CHK;

    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, 
            FN, FILE_MODE_READ, &reader);                                    CHK;

    ierr = VecLoad(vec, reader);                                             CHK;

    ierr = PetscViewerDestroy(&reader);                                      CHK;

    // get and print rhs information
    {
        PetscInt    n;
        ierr = VecGetSize(vec, &n);                                          CHK;
        ierr = PetscPrintf(PETSC_COMM_WORLD, 
                "Vector %s: size %d\n", name, n);                            CHK;
    }

    return 0;
}
