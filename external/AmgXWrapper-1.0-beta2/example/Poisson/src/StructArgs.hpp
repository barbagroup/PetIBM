/**
 * @file structArgs.cpp
 * @brief Definition of the structure StructArgs
 * @author Pi-Yueh Chuang (pychuang@gwu.edu)
 * @version alpha
 * @date 2015-11-06
 */
# include <petscsys.h>
# include <iostream>

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
        PetscInt            Nx,     // number of elements in x-direction
                            Ny,     // number of elements in y-direction
                            Nz = 0;     // number of elements in z-direction

        PetscInt            Nruns = 10; // number of runs of solving

        PetscBool           optFileBool,
                            VTKFileBool;

        char                mode[MAX_LEN],        // either AmgX, or PETSc
                            cfgFileName[MAX_LEN], // config file
                            optFileName[MAX_LEN], // output file
                            caseName[MAX_LEN],    // case name
                            VTKFileName[MAX_LEN]; // VTK output file name

        void print()
        {
            std::cout << "Case Name: " << caseName << std::endl;
            std::cout << "Nx: " << Nx << std::endl;
            std::cout << "Ny: " << Ny << std::endl;
            std::cout << "Nz: " << Nz << std::endl;
            std::cout << "Mode: " << mode << std::endl;
            std::cout << "Config File: " << cfgFileName << std::endl;
            std::cout << "Number of Solves: " << Nruns << std::endl;

            std::cout << "Output PETSc Log File ? " << optFileBool << std::endl;
            if (optFileBool == PETSC_TRUE)
                std::cout << "Output PETSc Log File Name: " 
                          << optFileName << std::endl;

            std::cout << "Output VTK file ? " << VTKFileBool << std::endl;
            if (VTKFileBool == PETSC_TRUE)
                std::cout << "Output VTK File Name: " 
                          << VTKFileName << std::endl;
        }

        PetscErrorCode getArgs()
        {
            PetscErrorCode      ierr;   // error codes returned by PETSc routines
            PetscBool           set;    // temporary booling variable

            ierr = PetscOptionsGetString(nullptr, nullptr, "-caseName", 
                    caseName, MAX_LEN, &set);                         CHKERRQ(ierr);
            CHKMSG(set, "caseName not yet set!");

            ierr = PetscOptionsGetString(nullptr, nullptr, "-mode", 
                    mode, MAX_LEN, &set);                             CHKERRQ(ierr);
            CHKMSG(set, "mode not yet set!");

            ierr = PetscOptionsGetString(nullptr, nullptr, "-cfgFileName", 
                    cfgFileName, MAX_LEN, &set);                      CHKERRQ(ierr);
            CHKMSG(set, "cfgFileName (configuration file) not yet set!");

            ierr = PetscOptionsGetString(nullptr, nullptr, "-optFileName", 
                    optFileName, MAX_LEN, &(optFileBool));       CHKERRQ(ierr);

            ierr = PetscOptionsGetString(nullptr, nullptr, "-VTKFileName", 
                    VTKFileName, MAX_LEN, &(VTKFileBool));       CHKERRQ(ierr);

            ierr = PetscOptionsGetInt(nullptr, nullptr, 
                    "-Nx", &Nx, &set);                                CHKERRQ(ierr);
            CHKMSG(set, "Nx not yet set!");

            ierr = PetscOptionsGetInt(nullptr, nullptr, 
                    "-Ny", &Ny, &set);                                CHKERRQ(ierr);
            CHKMSG(set, "Ny not yet set!");

            ierr = PetscOptionsGetInt(nullptr, nullptr, 
                    "-Nz", &Nz, &set);                                CHKERRQ(ierr);

            ierr = PetscOptionsGetInt(nullptr, nullptr, 
                    "-Nruns", &Nruns, &set);                          CHKERRQ(ierr);

            return 0;
        }
};
