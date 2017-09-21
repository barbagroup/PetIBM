/*
 * linsolver.cpp
 * Copyright (C) 2017 Pi-Yueh Chuang <pychuang@gwu.edu>
 *
 * Distributed under terms of the MIT license.
 */


// PetIBM
# include <petibm/linsolver.h>
# include <petibm/linsolverksp.h>

# ifdef HAVE_AMGX
# include <petibm/linsolveramgx.h>
# endif


namespace petibm
{
namespace linsolver
{

// factory function for creating a LinSolver
PetscErrorCode createLinSolver(
        const std::string &solverName, const std::string &configFile,
        const type::ExecuteType &type, type::LinSolver &solver)
{
    PetscFunctionBeginUser;

    switch (type)
    {
        case type::ExecuteType::CPU:
            solver = std::make_shared<LinSolverKSP>(solverName, configFile);
            break;
        case type::ExecuteType::GPU:
# ifdef HAVE_AMGX
            solver = std::make_shared<LinSolverAmgX>(solverName, configFile);
# else
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                    "AmgX solver is used, "
                    "while PetIBM is not compiled with AmgX.");
# endif
            break;
        default:
            break;
    }

    PetscFunctionReturn(0);
}

} // end of namespace linsolver
} // end of namespace petibm
