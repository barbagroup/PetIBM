/**
 * \file StructArgs.cpp
 * \brief definition of member functions in StructArgs. 
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \date 2017-06-01
 */


// class header
# include "StructArgs.hpp"


// macro for simplicity
# define CHKMSG(flag, message)                  \
    if (!flag) SETERRQ(PETSC_COMM_WORLD, 91, message);


// definition of print
PetscErrorCode StructArgs::print()
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    ierr = PetscPrintf(PETSC_COMM_WORLD, "Case Name: %s\n", caseName); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Nx: %d\n", Nx); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Ny: %d\n", Ny); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Nz: %d\n", Nz); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Mode: %s\n", mode); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Config File: %s\n", cfgFileName); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Number of Solves: %d\n", Nruns); CHKERRQ(ierr);

    ierr = PetscPrintf(PETSC_COMM_WORLD, "Output PETSc Log File ? %s\n", 
            optFileBool?"true":"false"); CHKERRQ(ierr);

    if (optFileBool) ierr = PetscPrintf(PETSC_COMM_WORLD,
            "Output PETSc Log File Name: %s\n", optFileName); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


// definition of checkHelp
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
        ierr = PetscPrintf(PETSC_COMM_WORLD,
                "\t-mode [PETSc or AmgX_GPU]\n"); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,
                "\t-cfgFileName [config file for solver]\n"); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,
                "\t-Nx [number of cells in x direction]\n"); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,
                "\t-Ny [number of cells in y direction]\n"); CHKERRQ(ierr);

        ierr = PetscPrintf(PETSC_COMM_WORLD, "Optional Parameters:\n"); CHKERRQ(ierr);

        ierr = PetscPrintf(PETSC_COMM_WORLD,
                "\t-Nz [number of cells in z direction for 3D cases]\n"); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,
                "\t-Nruns [number of runs in addition to warm-up run]\n"); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD, "\t-optFileName "
                "[file name for outputing PETSc performance log file]\n"); CHKERRQ(ierr);

        ierr = PetscFinalize(); CHKERRQ(ierr);

        exit(EXIT_SUCCESS);
    }

    PetscFunctionReturn(0);
}


// definition of getArgs
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

    ierr = PetscOptionsGetInt(nullptr, nullptr, "-Nx", &Nx, &set); CHKERRQ(ierr);
    CHKMSG(set, "Nx not yet set!");

    ierr = PetscOptionsGetInt(nullptr, nullptr, "-Ny", &Ny, &set); CHKERRQ(ierr);
    CHKMSG(set, "Ny not yet set!");

    ierr = PetscOptionsGetInt(nullptr, nullptr, "-Nz", &Nz, &set); CHKERRQ(ierr);

    ierr = PetscOptionsGetInt(nullptr, nullptr, "-Nruns", &Nruns, &set); CHKERRQ(ierr);

    ierr = PetscOptionsGetString(nullptr, nullptr, "-optFileName", 
            optFileName, PETSC_MAX_PATH_LEN, &optFileBool); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
