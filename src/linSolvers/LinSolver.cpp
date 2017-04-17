/*
 * LinSolver.cpp
 * Copyright (C) 2017 Pi-Yueh Chuang <pychuang@gwu.edu>
 *
 * Distributed under terms of the MIT license.
 */


// PetIBM
# include "LinSolver.h"
# include "KSPSolver.h"

# ifdef HAVE_AMGX
# include "AMGXSolver.h"
# endif


/** \copydoc createLinSolver. */
PetscErrorCode createLinSolver(const std::string &solverName, 
        const std::string &configFile, const types::ExecuteType &type,
        std::shared_ptr<LinSolver> &solver)
{
    PetscFunctionBeginUser;

    switch (type)
    {
        case types::ExecuteType::CPU:
            solver = std::shared_ptr<LinSolver>(
                    new KSPSolver(solverName, configFile));
            break;
        case types::ExecuteType::GPU:
# ifdef HAVE_AMGX
            solver = std::shared_ptr<LinSolver>(
                    new AMGXSolver(solverName, configFile));
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
