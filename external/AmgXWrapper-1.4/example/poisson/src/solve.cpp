/**
 * @file solve.cpp
 * @brief functions of solve.
 * @author Pi-Yueh Chuang (pychuang@gwu.edu)
 * @date 2016-02-04
 */


// PETSc
# include <petsctime.h>

// prototypes
# include "solve.hpp"


// definition of solve (KSP version)
PetscErrorCode solve(KSP &ksp, Mat &A, Vec &lhs, Vec &rhs, Vec &exact, Vec &err,
        StructArgs &args, PetscLogEvent &warmUpEvent, PetscLogEvent &solvingEvent)
{
    PetscFunctionBeginUser;

    PetscErrorCode          ierr;


    PetscLogDouble          totalTime = 0.,
                            avgTime;


    // tasks that will be carried out after each solving
    auto postSolving = [&ierr, &ksp, &lhs, &exact, &err](PetscLogDouble time)->PetscErrorCode
    {
        PetscFunctionBeginUser;

        KSPConvergedReason      reason; // to store the KSP convergence reason

        PetscInt                Niters; // iterations used to converge
        PetscInt                N; // Vector size

        PetscScalar             norm2,  // 2 norm of solution errors
                                normM;  // infinity norm of solution errors

        ierr = KSPGetConvergedReason(ksp, &reason); CHKERRQ(ierr);

        if (reason < 0) SETERRQ1(PETSC_COMM_WORLD,
                PETSC_ERR_CONV_FAILED, "Diverger reason: %d\n", reason);

        ierr = KSPGetIterationNumber(ksp, &Niters); CHKERRQ(ierr);

        // calculate norms of errors
        ierr = VecCopy(lhs, err); CHKERRQ(ierr);
        ierr = VecAXPY(err, -1.0, exact); CHKERRQ(ierr);
        ierr = VecNorm(err, NORM_2, &norm2); CHKERRQ(ierr);
        ierr = VecNorm(err, NORM_INFINITY, &normM); CHKERRQ(ierr);
        ierr = VecGetSize(err, &N); CHKERRQ(ierr);

        // print infromation
        ierr = PetscPrintf(PETSC_COMM_WORLD, "\tSolve Time: %f\n", time); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD, "\tL2-Norm: %g\n", (double)norm2/PetscSqrtReal(N)); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD, "\tMax-Norm: %g\n", (double)normM); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD, "\tIterations %D\n", Niters); CHKERRQ(ierr); 

        PetscFunctionReturn(0);
    };


    // tasks for soving with performance profiling
    auto solving = [&ierr, &A, &lhs, &rhs, &ksp, &postSolving](PetscLogEvent &event)->PetscLogDouble
    {
        PetscLogDouble          bgTime, 
                                edTime;

        // set all entries as zeros in the vector of unknows
        ierr = VecSet(lhs, 0.0); CHKERRQ(ierr);

        // begin logging this solving event
        ierr = PetscLogEventBegin(event, (PetscObject) A,
                (PetscObject) rhs, (PetscObject) lhs, 0); CHKERRQ(ierr);

        // obtain the begining time 
        ierr = PetscTime(&bgTime); CHKERRQ(ierr);

        // solve
        ierr = KSPSolve(ksp, rhs, lhs); CHKERRQ(ierr);

        // obtain the ending time 
        ierr = PetscTime(&edTime); CHKERRQ(ierr);

        // end logging this solving event
        ierr = PetscLogEventEnd(event, (PetscObject) A,
                (PetscObject) rhs, (PetscObject) lhs, 0); CHKERRQ(ierr);

        // post-solving tasks
        ierr = postSolving(edTime - bgTime); CHKERRQ(ierr);

        return edTime - bgTime;
    };


    // start a warm-up cycle
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Warm-up cycle ... \n"); CHKERRQ(ierr);
    solving(warmUpEvent);


    // start `Nruns` runs to get averaged solving time
    for(int i=0; i<args.Nruns; ++i)
    {
        ierr = PetscPrintf(PETSC_COMM_WORLD, "Run # %d ... \n", i); CHKERRQ(ierr);
        totalTime += solving(solvingEvent);
    }

    avgTime = totalTime / (double) args.Nruns;

    ierr = PetscPrintf(PETSC_COMM_WORLD,
            "\nThe averaged solving time is: %f\n", avgTime); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


// definition of solve (AmgX version)
PetscErrorCode solve(AmgXSolver &amgx, Mat &A, Vec &lhs, Vec &rhs, Vec &exact, Vec &err,
        StructArgs &args, PetscLogEvent &warmUpEvent, PetscLogEvent &solvingEvent)
{
    PetscFunctionBeginUser;

    PetscErrorCode          ierr;


    PetscLogDouble          totalTime = 0.,
                            avgTime;


    // collection of post-solving tasks
    auto postSolving = [&ierr, &amgx, &lhs, &exact, &err](PetscLogDouble time)->PetscErrorCode
    {
        PetscFunctionBeginUser;

        PetscInt                Niters; // iterations used to converge
        PetscInt                N; // Vector size

        PetscScalar             norm2,  // 2 norm of solution errors
                                normM;  // infinity norm of solution errors

        ierr = amgx.getIters(Niters); CHKERRQ(ierr);

        // calculate norms of errors
        ierr = VecCopy(lhs, err); CHKERRQ(ierr);
        ierr = VecAXPY(err, -1.0, exact); CHKERRQ(ierr);
        ierr = VecNorm(err, NORM_2, &norm2); CHKERRQ(ierr);
        ierr = VecNorm(err, NORM_INFINITY, &normM); CHKERRQ(ierr);
        ierr = VecGetSize(err, &N); CHKERRQ(ierr);

        // print infromation
        ierr = PetscPrintf(PETSC_COMM_WORLD, "\tSolve Time: %f\n", time); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD, "\tL2-Norm: %g\n", (double)norm2/PetscSqrtReal(N)); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD, "\tMax-Norm: %g\n", (double)normM); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD, "\tIterations %D\n", Niters); CHKERRQ(ierr); 

        PetscFunctionReturn(0);
    };


    // solving tasks with performance profiling
    auto solving = [&ierr, &A, &lhs, &rhs, &amgx, &postSolving](PetscLogEvent &event)->PetscLogDouble
    {
        PetscLogDouble          bgTime, 
                                edTime;

        // set all entries as zeros in the vector of unknows
        ierr = VecSet(lhs, 0.0); CHKERRQ(ierr);

        // beging logging this solving event
        ierr = PetscLogEventBegin(event, (PetscObject) A,
                (PetscObject) rhs, (PetscObject) lhs, 0); CHKERRQ(ierr);

        // obtain the begining time 
        ierr = PetscTime(&bgTime); CHKERRQ(ierr);

        // solve
        ierr = amgx.solve(lhs, rhs); CHKERRQ(ierr);

        // obtain the ending time 
        ierr = PetscTime(&edTime); CHKERRQ(ierr);

        // end logging this solving event
        ierr = PetscLogEventEnd(event, (PetscObject) A,
                (PetscObject) rhs, (PetscObject) lhs, 0); CHKERRQ(ierr);

        ierr = postSolving(edTime - bgTime); CHKERRQ(ierr);

        return edTime - bgTime;
    };


    // start a warm-up cycle
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Warm-up cycle ... \n"); CHKERRQ(ierr);
    solving(warmUpEvent);


    // start `Nruns` runs to get averaged solving time
    for(int i=0; i<args.Nruns; ++i)
    {
        ierr = PetscPrintf(PETSC_COMM_WORLD, "Run # %d ... \n", i); CHKERRQ(ierr);
        totalTime += solving(solvingEvent);
    }

    avgTime = totalTime / (double) args.Nruns;

    ierr = PetscPrintf(PETSC_COMM_WORLD,
            "\nThe averaged solving time is: %f\n", avgTime); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
