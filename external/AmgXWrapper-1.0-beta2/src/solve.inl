/**
 * @file solve.inl
 * @brief Definition of some member functions of the class AmgXSolver
 * @author Pi-Yueh Chuang (pychuang@gwu.edu)
 * @version alpha
 * @date 2016-01-08
 */
# include "AmgXSolver.hpp"


int AmgXSolver::solve(Vec &p, Vec &b)
{
    PetscErrorCode      ierr;

    if (redistScatter != nullptr)
    {
        ierr = VecScatterBegin(redistScatter, 
                b, redistRhs, INSERT_VALUES, SCATTER_FORWARD);               CHK;
        ierr = VecScatterEnd(redistScatter, 
                b, redistRhs, INSERT_VALUES, SCATTER_FORWARD);               CHK;

        ierr = VecScatterBegin(redistScatter, 
                p, redistLhs, INSERT_VALUES, SCATTER_FORWARD);               CHK;
        ierr = VecScatterEnd(redistScatter, 
                p, redistLhs, INSERT_VALUES, SCATTER_FORWARD);               CHK;

        if (gpuWorld != MPI_COMM_NULL)
        {
            ierr = solve_real(redistLhs, redistRhs);                         CHK;
        }
        MPI_Barrier(globalCpuWorld);

        ierr = VecScatterBegin(redistScatter, 
                redistLhs, p, INSERT_VALUES, SCATTER_REVERSE);               CHK;
        ierr = VecScatterEnd(redistScatter, 
                redistLhs, p, INSERT_VALUES, SCATTER_REVERSE);               CHK;
    }
    else
    {
        if (gpuWorld != MPI_COMM_NULL)
        {
            ierr = solve_real(p, b);                                         CHK;
        }
        MPI_Barrier(globalCpuWorld);
    }

    return 0;
}


/**
 * @brief Solve the linear problem based on given RHS and initial guess
 *
 * This function solve the linear problem. Users need to set matrix A before
 * calling this function. The result will be saved back to the PETSc Vec 
 * containing initial guess.
 *
 * It lacks mechanisms to check whether necessary initialization and setting
 * matrix A are done first.
 *
 * @param p Unknowns vector. The values passed in will be an initial guess.
 * @param b Right-hand-side vector.
 *
 * @return Currently meaningless. May be error codes in the future.
 */
int AmgXSolver::solve_real(Vec &p, Vec &b)
{
    PetscErrorCode      ierr;
    double             *unks,
                       *rhs;
    int                 size;
    AMGX_SOLVE_STATUS   status;

    // get size of local vector (p and b should have the same local size)
    ierr = VecGetLocalSize(p, &size);                                        CHK;

    // get pointers to the raw data of local vectors
    ierr = VecGetArray(p, &unks);                                            CHK;
    ierr = VecGetArray(b, &rhs);                                             CHK;

    // upload vectors to AmgX
    AMGX_vector_upload(AmgXP, size, 1, unks);
    AMGX_vector_upload(AmgXRHS, size, 1, rhs);

    // solve
    MPI_Barrier(gpuWorld);
    AMGX_solver_solve(solver, AmgXRHS, AmgXP);

    // get the status of the solver
    AMGX_solver_get_status(solver, &status);

    // check whether the solver successfully solve the problem
    if (status != AMGX_SOLVE_SUCCESS)
        std::cout << "AmgX solver failed to solve the problem! "
                  << "Solver status: " << status << std::endl;

    // download data from device
    AMGX_vector_download(AmgXP, unks);

    // restore PETSc vectors
    ierr = VecRestoreArray(p, &unks);                                        CHK;
    ierr = VecRestoreArray(b, &rhs);                                         CHK;

    return 0;
}

