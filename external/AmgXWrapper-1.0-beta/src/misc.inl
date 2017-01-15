/**
 * @file misc.cpp
 * @brief Definition of some member functions of the class AmgXSolver
 * @author Pi-Yueh Chuang (pychuang@gwu.edu)
 * @version alpha
 * @date 2016-01-08
 */
# include "AmgXSolver.hpp"


/**
 * @brief Setting up AmgX mode.
 *
 * Convert a STL string to AmgX mode and then store this mode. For the usage of 
 * AmgX modes, please refer to AmgX manual.
 *
 * @param _mode A STL string describing the mode.
 *
 * @return Currently meaningless. May be error codes in the future.
 */
int AmgXSolver::setMode(const std::string &_mode)
{
    if (_mode == "hDDI")
        mode = AMGX_mode_hDDI;
    else if (_mode == "hDFI")
        mode = AMGX_mode_hDFI;
    else if (_mode == "hFFI")
        mode = AMGX_mode_hFFI;
    else if (_mode == "dDDI")
        mode = AMGX_mode_dDDI;
    else if (_mode == "dDFI")
        mode = AMGX_mode_dDFI;
    else if (_mode == "dFFI")
        mode = AMGX_mode_dFFI;
    else
    {
        std::cerr << "error: " 
                  << _mode << " is not an available mode." << std::endl;
        exit(0);
    }
    return 0;
}


/**
 * @brief Get the number of iterations of last solve phase
 *
 * @return number of iterations
 */
int AmgXSolver::getIters()
{
    int         iter = 0;

    if (gpuProc == 0)
        AMGX_solver_get_iterations_number(solver, &iter); 

    return iter;
}


/**
 * @brief Get the residual at a specific iteration in last solve phase
 *
 * @param iter a specific iteration during the last solve phase
 *
 * @return residual
 */
double AmgXSolver::getResidual(const int &iter)
{
    double      res = 0;

    if (gpuProc == 0)
        AMGX_solver_get_iteration_residual(solver, iter, 0, &res);

    return res;
}


/**
 * @brief A printing function, using stdout, needed by AmgX.
 *
 * @param msg C-style string
 * @param length The length of the string
 */
void AmgXSolver::print_callback(const char *msg, int length)
{
    std::cout << msg;
}


/**
 * @brief A printing function, print nothing, needed by AmgX.
 *
 * @param msg C-style string
 * @param length The length of the string
 */
void AmgXSolver::print_none(const char *msg, int length) { }


/**
 * @brief Get the memory usage on devices
 *
 * @return error codes in the future. 
 */
int AmgXSolver::getMemUsage()
{
    size_t free_byte,
           total_byte;

    if (gpuProc == 0)
    {
        CHECK(cudaMemGetInfo(&free_byte, &total_byte));

        std::cout << "myGlobalRank: " << myGlobalRank << " "
                  << free_byte / 1024.0 / 1024.0 << " MB "
                  << " / " << total_byte / 1024.0 / 1024.0 << " MB " << std::endl;
    }

    return 0;
}


