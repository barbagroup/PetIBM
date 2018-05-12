/**
 * @file structArgs.cpp
 * @brief definition of the member functions of structure StructArgs
 * @author Pi-Yueh Chuang (pychuang@gwu.edu)
 * @version alpha
 * @date 2015-11-06
 */


// class definition
# include "StructArgs.hpp"

// macro for simplicity
# define CHKMSG(flag, message)                  \
    if (!flag) SETERRQ(PETSC_COMM_WORLD, 91, message);


// definition of StructArgs::print
PetscErrorCode StructArgs::print()
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    ierr = PetscPrintf(PETSC_COMM_WORLD, "Case Name: %s\n", caseName); CHKERRQ(ierr);

    ierr = PetscPrintf(PETSC_COMM_WORLD, "Mode: %s\n", mode); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Config File: %s\n", cfgFileName); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Matrix File: %s\n", matrixFileName); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "RHS File: %s\n", rhsFileName); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Exact Solution File: %s\n", exactFileName); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Number of Solves: %d\n", Nruns); CHKERRQ(ierr);

    if (optFileBool == PETSC_TRUE) ierr = PetscPrintf(PETSC_COMM_WORLD,
            "Output PETSc Log File Name: %s\n", optFileName); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


// definition of StructArgs::checkHelp
PetscErrorCode StructArgs::checkHelp()
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;
    PetscBool           help = PETSC_FALSE;

    ierr = PetscOptionsGetBool(nullptr, nullptr, "-print", &help, nullptr); CHKERRQ(ierr);

    if (help)
    {
        ierr = PetscPrintf(PETSC_COMM_WORLD, "Necessary Parameters:\n"); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD, "\t-caseName [string]\n"); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD, "\t-mode [PETSc, AmgX_CPU, AmgX_GPU]\n"); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD, "\t-cfgFileName [config file for solver]\n"); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD, "\t-matrixFileName [file for matrix]\n"); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD, "\t-rhsFileName [file for RHS vector]\n"); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,
                "\t-exactFileName [file for exact solution vector]\n"); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD, "\n"); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD, "Optional Parameters:\n"); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD, 
                "\t-optFileName [file name for outputing "
                "PETSc performance log file]\n"); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD, "\t-Nruns "
                "[number of runs in addition to warm-up run]\n"); CHKERRQ(ierr);

        ierr = PetscFinalize(); CHKERRQ(ierr);

        exit(EXIT_SUCCESS);
    }

    PetscFunctionReturn(0);
}


// definition of StructArgs::getArgs
PetscErrorCode StructArgs::getArgs()
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;   // error codes returned by PETSc routines
    PetscBool           set;    // temporary booling variable

    ierr = checkHelp(); CHKERRQ(ierr);

    ierr = PetscOptionsGetString(nullptr, nullptr, "-caseName", 
            caseName, PETSC_MAX_PATH_LEN, &set); CHKERRQ(ierr);
    CHKMSG(set, "caseName not yet set!");

    ierr = PetscOptionsGetString(nullptr, nullptr, "-mode", 
            mode, PETSC_MAX_PATH_LEN, &set); CHKERRQ(ierr);
    CHKMSG(set, "mode not yet set!");

    ierr = PetscOptionsGetString(nullptr, nullptr, "-cfgFileName", 
            cfgFileName, PETSC_MAX_PATH_LEN, &set); CHKERRQ(ierr);
    CHKMSG(set, "cfgFileName (configuration file) not yet set!");

    ierr = PetscOptionsGetString(nullptr, nullptr, "-optFileName", 
            optFileName, PETSC_MAX_PATH_LEN, &optFileBool); CHKERRQ(ierr);

    ierr = PetscOptionsGetString(nullptr, nullptr, "-matrixFileName", 
            matrixFileName, PETSC_MAX_PATH_LEN, &set); CHKERRQ(ierr);
    CHKMSG(set, "matrixFileName (binary matrix file) not yet set!");

    ierr = PetscOptionsGetString(nullptr, nullptr, "-rhsFileName", 
            rhsFileName, PETSC_MAX_PATH_LEN, &set); CHKERRQ(ierr);
    CHKMSG(set, "rhsFileName (binary matrix file) not yet set!");

    ierr = PetscOptionsGetString(nullptr, nullptr, "-exactFileName", 
            exactFileName, PETSC_MAX_PATH_LEN, &set); CHKERRQ(ierr);
    CHKMSG(set, "exactFileName (binary matrix file) not yet set!");

    ierr = PetscOptionsGetInt(nullptr, nullptr, "-Nruns", &Nruns, &set); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
