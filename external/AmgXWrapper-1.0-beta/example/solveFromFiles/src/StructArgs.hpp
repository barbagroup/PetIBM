/**
 * @file structArgs.cpp
 * @brief Definition of the structure StructArgs
 * @author Pi-Yueh Chuang (pychuang@gwu.edu)
 * @version alpha
 * @date 2015-11-06
 */
# pragma once

# define MAX_LEN PETSC_MAX_PATH_LEN

# define CHKMSG(flag, message)                  \
    if (!flag)                                  \
    {                                           \
        PetscPrintf(PETSC_COMM_WORLD, message); \
        PetscPrintf(PETSC_COMM_WORLD, "\n");    \
        exit(EXIT_FAILURE);                     \
    }


class StructArgs
{
    public: 
        PetscInt            Nruns = 10; // number of runs of solving

        PetscBool           optFileBool;

        char                mode[MAX_LEN],        // either AmgX, or PETSc
                            cfgFileName[MAX_LEN], // config file
                            optFileName[MAX_LEN], // output file
                            matrixFileName[MAX_LEN], // output file
                            rhsFileName[MAX_LEN], // output file
                            exactFileName[MAX_LEN],  // output file
                            caseName[MAX_LEN];    // case name

        void print()
        {
            std::cout << "Case Name: " << caseName << std::endl;
            std::cout << "Mode: " << mode << std::endl;
            std::cout << "Config File: " << cfgFileName << std::endl;
            std::cout << "Matrix File: " << matrixFileName << std::endl;
            std::cout << "RHS File: " << rhsFileName << std::endl;
            std::cout << "Exact Solution File: " << exactFileName << std::endl;
            std::cout << "Number of Solves: " << Nruns << std::endl;

            std::cout << "Output PETSc Log File ? " << optFileBool << std::endl;
            if (optFileBool == PETSC_TRUE)
                std::cout << "Output PETSc Log File Name: " 
                          << optFileName << std::endl;
        }


        PetscErrorCode checkHelp()
        {
            PetscErrorCode      ierr;
            PetscBool           help,
                                set;

            ierr = PetscOptionsGetBool(nullptr, nullptr, "-print", &help, &set); CHK;

            if (set == PETSC_TRUE)
            {
                ierr = PetscPrintf(PETSC_COMM_WORLD, 
                        "Necessary Parameters:\n");                          CHK;
                ierr = PetscPrintf(PETSC_COMM_WORLD, 
                        "\t-caseName [string]\n");                           CHK;
                ierr = PetscPrintf(PETSC_COMM_WORLD, 
                        "\t-mode [PETSc or AmgX]\n");                        CHK;
                ierr = PetscPrintf(PETSC_COMM_WORLD, 
                        "\t-cfgFileName [config file for solver]\n");        CHK;
                ierr = PetscPrintf(PETSC_COMM_WORLD, 
                        "\t-matrixFileName [file for matrix]\n");            CHK;
                ierr = PetscPrintf(PETSC_COMM_WORLD, 
                        "\t-rhsFileName [file for RHS vector]\n");           CHK;
                ierr = PetscPrintf(PETSC_COMM_WORLD, 
                        "\t-exactFileName [file for exact solution vector]\n"); CHK;

                ierr = PetscPrintf(PETSC_COMM_WORLD, "\n");                  CHK;

                ierr = PetscPrintf(PETSC_COMM_WORLD, 
                        "Optional Parameters:\n");                           CHK;
                ierr = PetscPrintf(PETSC_COMM_WORLD, 
                        "\t-optFileName [file name for outputing "
                        "PETSc performance log file]\n");                    CHK;
                ierr = PetscPrintf(PETSC_COMM_WORLD, 
                        "\t-Nruns "
                        "[number of runs in addition to warm-up run]\n");    CHK;

                exit(EXIT_SUCCESS);
            }

            return 0;
        }

        PetscErrorCode getArgs()
        {
            PetscErrorCode      ierr;   // error codes returned by PETSc routines
            PetscBool           set;    // temporary booling variable

            ierr = checkHelp();                                              CHK;

            ierr = PetscOptionsGetString(nullptr, nullptr, "-caseName", 
                    caseName, MAX_LEN, &set);                                CHK;
            CHKMSG(set, "caseName not yet set!");

            ierr = PetscOptionsGetString(nullptr, nullptr, "-mode", 
                    mode, MAX_LEN, &set);                                    CHK;
            CHKMSG(set, "mode not yet set!");

            ierr = PetscOptionsGetString(nullptr, nullptr, "-cfgFileName", 
                    cfgFileName, MAX_LEN, &set);                             CHK;
            CHKMSG(set, "cfgFileName (configuration file) not yet set!");

            ierr = PetscOptionsGetString(nullptr, nullptr, "-optFileName", 
                    optFileName, MAX_LEN, &optFileBool);                     CHK;

            ierr = PetscOptionsGetString(nullptr, nullptr, "-matrixFileName", 
                    matrixFileName, MAX_LEN, &set);                          CHK;
            CHKMSG(set, "matrixFileName (binary matrix file) not yet set!");

            ierr = PetscOptionsGetString(nullptr, nullptr, "-rhsFileName", 
                    rhsFileName, MAX_LEN, &set);                             CHK;
            CHKMSG(set, "rhsFileName (binary matrix file) not yet set!");

            ierr = PetscOptionsGetString(nullptr, nullptr, "-exactFileName", 
                    exactFileName, MAX_LEN, &set);                           CHK;
            CHKMSG(set, "exactFileName (binary matrix file) not yet set!");

            ierr = PetscOptionsGetInt(nullptr, nullptr, 
                    "-Nruns", &Nruns, &set);                                 CHK;

            return 0;
        }
};
