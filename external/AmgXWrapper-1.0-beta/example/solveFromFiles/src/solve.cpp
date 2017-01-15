/**
 * @file solve.cpp
 * @brief functions of solving behavior
 * @author Pi-Yueh Chuang (pychuang@gwu.edu)
 * @version alpha
 * @date 2016-02-04
 */
# include "headers.hpp"


PetscErrorCode solve(KSP &ksp, Mat &A, Vec &u, Vec &rhs, Vec &u_exact, Vec &err,
        StructArgs &args, PetscLogEvent &warmUpEvent, PetscLogEvent &solvingEvent)
{

    PetscErrorCode          ierr;


    PetscLogDouble          totalTime = 0.,
                            avgTime;


    auto                    postSolving =
        [&ierr, &ksp, &u, &u_exact, &err] (PetscLogDouble time) -> PetscErrorCode
        {
            KSPConvergedReason      reason; // to store the KSP convergence reason

            PetscInt                Niters; // iterations used to converge

            PetscScalar               norm2,  // L2 norm of solution errors
                                    normM;  // infinity norm of solution errors

            ierr = KSPGetConvergedReason(ksp, &reason);                      CHK;

            if (reason < 0)
            {
                ierr = PetscPrintf(PETSC_COMM_WORLD, 
                        "\nDiverged: %d\n", reason);                         CHK;
                exit(EXIT_FAILURE);
            }

            ierr = KSPGetIterationNumber(ksp, &Niters);                      CHK;

            // calculate norms of errors
            ierr = VecCopy(u, err);                                          CHK;
            ierr = VecAXPY(err, -1.0, u_exact);                              CHK;
            ierr = VecNorm(err, NORM_2, &norm2);                             CHK;
            ierr = VecNorm(err, NORM_INFINITY, &normM);                      CHK;

            // print infromation
            ierr = PetscPrintf(PETSC_COMM_WORLD, 
                    "\tSolve Time: %f\n", time);                  CHK;
            ierr = PetscPrintf(PETSC_COMM_WORLD, 
                    "\tL2-Norm: %g\n", (double)norm2);                       CHK;
            ierr = PetscPrintf(PETSC_COMM_WORLD, 
                    "\tMax-Norm: %g\n", (double)normM);                      CHK;
            ierr = PetscPrintf(PETSC_COMM_WORLD, 
                    "\tIterations %D\n", Niters);                            CHK; 

            return 0;
        };


    auto                    solving = 
        [&ierr, &A, &u, &rhs, &ksp, &postSolving] (PetscLogEvent &event) -> PetscLogDouble
        {
            PetscLogDouble          bgTime, 
                                    edTime;

            // set all entries as zeros in the vector of unknows
            ierr = MPI_Barrier(PETSC_COMM_WORLD);                            CHK;
            ierr = VecSet(u, 0.0);                                           CHK;

            // beging logging this solving event
            ierr = PetscLogEventBegin(event, 
                (PetscObject) A, (PetscObject) rhs, (PetscObject) u, 0);     CHK;

            // obtain the begining time 
            ierr = PetscTime(&bgTime);                                       CHK;

            // solve
            ierr = KSPSolve(ksp, rhs, u);                                    CHK;

            // obtain the ending time 
            ierr = PetscTime(&edTime);                                       CHK;

            // end logging this solving event
            ierr = PetscLogEventEnd(event, 
                (PetscObject) A, (PetscObject) rhs, (PetscObject) u, 0);     CHK;

            ierr = postSolving(edTime - bgTime);                             CHK;

            return edTime - bgTime;
        };




    // start a warm-up cycle
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Warm-up cycle ... \n");            CHK;
    solving(warmUpEvent);


    // start `Nruns` runs to get averaged solving time
    for(int i=0; i<args.Nruns; ++i)
    {
        ierr = PetscPrintf(PETSC_COMM_WORLD, "Run # %d ... \n", i);          CHK;
        totalTime += solving(solvingEvent);
    }

    avgTime = totalTime / (double) args.Nruns;

    ierr = PetscPrintf(PETSC_COMM_WORLD,
            "\nThe averaged solving time is: %f\n", avgTime);                  CHK;

    // destroy KSP
    ierr = KSPDestroy(&ksp);                                                 CHK;


    return 0;
}


PetscErrorCode solve(AmgXSolver &amgx, Mat &A, Vec &u, Vec &rhs, Vec &u_exact, Vec &err,
        StructArgs &args, PetscLogEvent &warmUpEvent, PetscLogEvent &solvingEvent)
{

    PetscErrorCode          ierr;


    PetscLogDouble          totalTime = 0.,
                            avgTime;


    auto                    postSolving =
        [&ierr, &amgx, &u, &u_exact, &err] (PetscLogDouble time) -> PetscErrorCode
        {
            PetscInt                Niters; // iterations used to converge

            PetscScalar               norm2,  // L2 norm of solution errors
                                    normM;  // infinity norm of solution errors

            Niters = amgx.getIters();

            // calculate norms of errors
            ierr = VecCopy(u, err);                                          CHK;
            ierr = VecAXPY(err, -1.0, u_exact);                              CHK;
            ierr = VecNorm(err, NORM_2, &norm2);                             CHK;
            ierr = VecNorm(err, NORM_INFINITY, &normM);                      CHK;

            // print infromation
            ierr = PetscPrintf(PETSC_COMM_WORLD, 
                    "\tSolve Time: %f\n", time);                  CHK;
            ierr = PetscPrintf(PETSC_COMM_WORLD, 
                    "\tL2-Norm: %g\n", (double)norm2);                       CHK;
            ierr = PetscPrintf(PETSC_COMM_WORLD, 
                    "\tMax-Norm: %g\n", (double)normM);                      CHK;
            ierr = PetscPrintf(PETSC_COMM_WORLD, 
                    "\tIterations %D\n", Niters);                            CHK; 

            return 0;
        };


    auto                    solving = 
        [&ierr, &A, &u, &rhs, &amgx, &postSolving] (PetscLogEvent &event) -> PetscLogDouble
        {
            PetscLogDouble          bgTime, 
                                    edTime;

            // set all entries as zeros in the vector of unknows
            ierr = MPI_Barrier(PETSC_COMM_WORLD);                            CHK;
            ierr = VecSet(u, 0.0);                                           CHK;

            // beging logging this solving event
            ierr = PetscLogEventBegin(event, 
                (PetscObject) A, (PetscObject) rhs, (PetscObject) u, 0);     CHK;

            // obtain the begining time 
            ierr = PetscTime(&bgTime);                                       CHK;

            // solve
            ierr = amgx.solve(u, rhs);                                       CHK;

            // obtain the ending time 
            ierr = PetscTime(&edTime);                                       CHK;

            // end logging this solving event
            ierr = PetscLogEventEnd(event, 
                (PetscObject) A, (PetscObject) rhs, (PetscObject) u, 0);     CHK;

            ierr = postSolving(edTime - bgTime);                             CHK;

            return edTime - bgTime;
        };


    // start a warm-up cycle
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Warm-up cycle ... \n");            CHK;
    solving(warmUpEvent);


    // start `Nruns` runs to get averaged solving time
    for(int i=0; i<args.Nruns; ++i)
    {
        ierr = PetscPrintf(PETSC_COMM_WORLD, "Run # %d ... \n", i);          CHK;
        totalTime += solving(solvingEvent);
    }

    avgTime = totalTime / (double) args.Nruns;

    ierr = PetscPrintf(PETSC_COMM_WORLD,
            "\nThe averaged solving time is: %f\n", avgTime);                  CHK;


    return 0;
}
