/**
 * @file StructArgs.hpp
 * @brief definition of the structure StructArgs
 * @author Pi-Yueh Chuang (pychuang@gwu.edu)
 * @date 2015-11-06
 */


# pragma once

// PETSc
# include <petscsys.h>


/** \brief struct for command-line arguments. */
class StructArgs
{
public: 

    PetscInt            Nx,         // number of elements in x-direction
                        Ny,         // number of elements in y-direction
                        Nz = 0;     // number of elements in z-direction

    PetscInt            Nruns = 10; // number of runs of solving

    PetscBool           optFileBool; // indicates if we will output a performance log.

    char                mode[PETSC_MAX_PATH_LEN],        // either AmgX_GPU, or PETSc
                        cfgFileName[PETSC_MAX_PATH_LEN], // config file
                        optFileName[PETSC_MAX_PATH_LEN], // output file
                        caseName[PETSC_MAX_PATH_LEN],    // case name
                        VTKFileName[PETSC_MAX_PATH_LEN]; // VTK output file name


    /** \brief default constructor. */
    StructArgs() = default;

    /** \brief default destructor. */
    ~StructArgs() = default;

    /**
     * \brief print information based on the input from command line.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode print();

    /**
     * \brief check if users have `-help` in command-line argument.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode checkHelp();

    /**
     * \brief get command-line arguments.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode getArgs();
};
