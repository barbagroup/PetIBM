/**
 * @file AmgXSolver.cpp
 * @brief Definition of member functions of the class AmgXSolver
 * @author Pi-Yueh Chuang (pychuang@gwu.edu)
 * @version alpha
 * @date 2015-09-01
 */
# include "AmgXSolver.hpp"


/**
 * @brief Initialization of the static member -- count
 *
 * count is used to count the number of instances. The fisrt instance is 
 * responsable for initializing AmgX library and the resource instance.
 */
int AmgXSolver::count = 0;


/**
 * @brief Initialization of the static member -- rsrc
 *
 * Due to the design of AmgX library, using more than one resource instance may
 * cause some problems. So make the resource instance as a static member to keep
 * only one instance there.
 */
AMGX_resources_handle AmgXSolver::rsrc = nullptr;

# include "init.inl"
# include "setA.inl"
# include "solve.inl"
# include "misc.inl"

