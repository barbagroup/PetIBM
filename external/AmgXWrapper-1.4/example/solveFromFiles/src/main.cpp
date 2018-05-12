/**
 * @file main.cpp
 * @brief An example and benchmark of AmgX and PETSc
 * @author Pi-Yueh Chuang (pychuang@gwu.edu)
 * @version beta
 * @date 2015-02-01
 */


// STL
# include <cstring>

// PETSc
# include <petscsys.h>
# include <petscmat.h>
# include <petscvec.h>
# include <petscksp.h>

// AmgXWrapper
# include <AmgXSolver.hpp>

// headers
# include "StructArgs.hpp"
# include "io.hpp"
# include "createKSP.hpp"
# include "solve.hpp"


int main(int argc, char **argv)
{
    StructArgs          args;   // a structure containing CMD arguments


    Vec                 lhs,    // unknowns
                        rhs,    // RHS
                        exact,  // exact solution
                        err;    // errors

    Mat                 A;      // coefficient matrix

    KSP                 ksp;    // PETSc KSP solver instance

    AmgXSolver          amgx;   // AmgX wrapper instance

    PetscErrorCode      ierr;   // error codes returned by PETSc routines

    PetscMPIInt         size,   // MPI size
                        myRank; // rank of current process

    PetscClassId        solvingID,
                        warmUpID;

    PetscLogEvent       solvingEvent,
                        warmUpEvent;





    // initialize PETSc and MPI
    ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);
    ierr = PetscLogDefaultBegin(); CHKERRQ(ierr);


    // obtain the rank and size of MPI
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &myRank); CHKERRQ(ierr);


    // get parameters from command-line arguments
    ierr = args.getArgs(); CHKERRQ(ierr);

    // print case information
    ierr = printHeader(args); CHKERRQ(ierr);

    // create matrix A and load from file
    ierr = readMat(A, args.matrixFileName, "A"); CHKERRQ(ierr);

    // create vector rhs and load from file
    ierr = readVec(rhs, args.rhsFileName, "RHS"); CHKERRQ(ierr);

    // create vector exact and load from file
    ierr = readVec(exact, args.exactFileName, "exact solution"); CHKERRQ(ierr);


    // create vectors based on the matrix and set to zeros
    {
        ierr = MatCreateVecs(A, &lhs, nullptr); CHKERRQ(ierr);
        ierr = PetscObjectSetName((PetscObject) lhs, "unknowns"); CHKERRQ(ierr);
        ierr = VecSet(lhs, 0.0); CHKERRQ(ierr);
    }

    // create vectors err and set to zeros
    {
        ierr = VecDuplicate(lhs, &err); CHKERRQ(ierr);
        ierr = PetscObjectSetName((PetscObject) err, "errors"); CHKERRQ(ierr);
        ierr = VecSet(err, 0.0); CHKERRQ(ierr);
    }


    // register a PETSc event for warm-up and solving
    ierr = PetscClassIdRegister("SolvingClass", &solvingID); CHKERRQ(ierr);
    ierr = PetscClassIdRegister("WarmUpClass", &warmUpID); CHKERRQ(ierr);
    ierr = PetscLogEventRegister("Solving", solvingID, &solvingEvent); CHKERRQ(ierr);
    ierr = PetscLogEventRegister("WarmUp", warmUpID, &warmUpEvent); CHKERRQ(ierr);

    // create a solver and solve based whether it is AmgX or PETSc
    if (std::strcmp(args.mode, "PETSc") == 0) // PETSc mode
    {
        ierr = createKSP(ksp, A, args.cfgFileName); CHKERRQ(ierr);

        ierr = solve(ksp, A, lhs, rhs, exact, err, 
                args, warmUpEvent, solvingEvent); CHKERRQ(ierr);

        // destroy KSP
        ierr = KSPDestroy(&ksp); CHKERRQ(ierr);

    }
    else // AmgX mode
    {
        // AmgX GPU mode
        if (std::strcmp(args.mode, "AmgX_GPU") == 0)
            amgx.initialize(PETSC_COMM_WORLD, "dDDI", args.cfgFileName);
        // AmgX CPU mode (not supported yet)
        //else if (std::strcmp(args.mode, "AmgX_CPU") == 0)
        //    amgx.initialize(PETSC_COMM_WORLD, "hDDI", args.cfgFileName);
        else SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_ARG_UNKNOWN_TYPE,
                    "Invalid mode: %s\n", args.mode); CHKERRQ(ierr);


        ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);
        ierr = amgx.setA(A); CHKERRQ(ierr);

        ierr = solve(amgx, A, lhs, rhs, exact, err, 
                args, warmUpEvent, solvingEvent); CHKERRQ(ierr);

        // destroy solver
        ierr = amgx.finalize(); CHKERRQ(ierr);

    }


    // output a file for petsc performance
    if (args.optFileBool == PETSC_TRUE)
    {
        PetscViewer         viewer; // PETSc viewer

        std::strcat(args.optFileName ,".log");

        ierr = PetscViewerASCIIOpen(
                PETSC_COMM_WORLD, args.optFileName, &viewer); CHKERRQ(ierr);
        ierr = PetscLogView(viewer); CHKERRQ(ierr);
        ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
    }
    

    // destroy vectors, matrix, dmda
    {
        ierr = VecDestroy(&lhs); CHKERRQ(ierr);
        ierr = VecDestroy(&rhs); CHKERRQ(ierr);
        ierr = VecDestroy(&exact); CHKERRQ(ierr);
        ierr = VecDestroy(&err); CHKERRQ(ierr);
        ierr = MatDestroy(&A); CHKERRQ(ierr);
    }


    // printer a footer
    ierr = printFooter(args); CHKERRQ(ierr);

    // finalize PETSc
    ierr = PetscFinalize(); CHKERRQ(ierr);

    return 0;
}
