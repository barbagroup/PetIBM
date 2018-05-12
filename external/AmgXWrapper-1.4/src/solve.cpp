/**
 * @file solve.cpp
 * @brief definition of member functions regarding to solving in AmgXSolver.
 * @author Pi-Yueh Chuang (pychuang@gwu.edu)
 * @date 2016-01-08
 */


// AmgXWrapper
# include "AmgXSolver.hpp"


// definition of AmgXSolver::solve
PetscErrorCode AmgXSolver::solve(Vec &p, Vec &b)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    if (globalSize != gpuWorldSize)
    {
        ierr = VecScatterBegin(scatterRhs, 
                b, redistRhs, INSERT_VALUES, SCATTER_FORWARD); CHK;
        ierr = VecScatterBegin(scatterLhs, 
                p, redistLhs, INSERT_VALUES, SCATTER_FORWARD); CHK;

        ierr = VecScatterEnd(scatterRhs, 
                b, redistRhs, INSERT_VALUES, SCATTER_FORWARD); CHK;
        ierr = VecScatterEnd(scatterLhs, 
                p, redistLhs, INSERT_VALUES, SCATTER_FORWARD); CHK;

        if (gpuWorld != MPI_COMM_NULL)
        {
            ierr = solve_real(redistLhs, redistRhs); CHK;
        }
        ierr = MPI_Barrier(globalCpuWorld); CHK;

        ierr = VecScatterBegin(scatterLhs, 
                redistLhs, p, INSERT_VALUES, SCATTER_REVERSE); CHK;
        ierr = VecScatterEnd(scatterLhs, 
                redistLhs, p, INSERT_VALUES, SCATTER_REVERSE); CHK;
    }
    else
    {
        if (gpuWorld != MPI_COMM_NULL)
        {
            ierr = solve_real(p, b); CHK;
        }
        ierr = MPI_Barrier(globalCpuWorld); CHK;
    }

    PetscFunctionReturn(0);
}


// definition of AmgXSolver::solve_real
PetscErrorCode AmgXSolver::solve_real(Vec &p, Vec &b)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    double              *unks,
                        *rhs;

    int                 size;

    AMGX_SOLVE_STATUS   status;

    // get size of local vector (p and b should have the same local size)
    ierr = VecGetLocalSize(p, &size); CHK;

    // get pointers to the raw data of local vectors
    ierr = VecGetArray(p, &unks); CHK;
    ierr = VecGetArray(b, &rhs); CHK;

    // upload vectors to AmgX
    AMGX_vector_upload(AmgXP, size, 1, unks);
    AMGX_vector_upload(AmgXRHS, size, 1, rhs);

    // solve
    ierr = MPI_Barrier(gpuWorld); CHK;
    AMGX_solver_solve(solver, AmgXRHS, AmgXP);

    // get the status of the solver
    AMGX_solver_get_status(solver, &status);

    // check whether the solver successfully solve the problem
    if (status != AMGX_SOLVE_SUCCESS) SETERRQ1(globalCpuWorld, 
            PETSC_ERR_CONV_FAILED, "AmgX solver failed to solve the system! "
            "The error code is %d.\n", status);

    // download data from device
    AMGX_vector_download(AmgXP, unks);

    // restore PETSc vectors
    ierr = VecRestoreArray(p, &unks); CHK;
    ierr = VecRestoreArray(b, &rhs); CHK;

    PetscFunctionReturn(0);
}
