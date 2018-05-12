/**
 * @file StructArgs.hpp
 * @brief An example and benchmark of AmgX and PETSc
 * @author Pi-Yueh Chuang (pychuang@gwu.edu)
 * @version beta
 * @date 2015-02-01
 */


# pragma once

// PETSc
# include <petscsys.h>


/** \brief a struct for command-line arguments. */
class StructArgs
{
public: 

    PetscInt            Nruns = 10; // number of runs of solving

    PetscBool           optFileBool;

    char                mode[PETSC_MAX_PATH_LEN],        // either AmgX, or PETSc
                        cfgFileName[PETSC_MAX_PATH_LEN], // config file
                        optFileName[PETSC_MAX_PATH_LEN], // output file
                        matrixFileName[PETSC_MAX_PATH_LEN], // output file
                        rhsFileName[PETSC_MAX_PATH_LEN], // output file
                        exactFileName[PETSC_MAX_PATH_LEN],  // output file
                        caseName[PETSC_MAX_PATH_LEN];    // case name

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
     * \brief check if users has `-help` in command-line argument.
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
