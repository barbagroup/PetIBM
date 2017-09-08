/*
 * linsolver.cpp
 * Copyright (C) 2017 Pi-Yueh Chuang <pychuang@gwu.edu>
 *
 * Distributed under terms of the MIT license.
 */


// PetIBM
# include "petibm/linsolver.h"
# include "petibm/kspsolver.h"

# ifdef HAVE_AMGX
# include "petibm/amgxsolver.h"
# endif


namespace petibm
{
namespace linsolvers
{

/** \copydoc createLinSolver. */
PetscErrorCode createLinSolver(const std::string &solverName, 
        const std::string &configFile,
        const petibm::utilities::types::ExecuteType &type,
        std::shared_ptr<LinSolver> &solver)
{
    PetscFunctionBeginUser;

    switch (type)
    {
        case petibm::utilities::types::ExecuteType::CPU:
            solver = std::shared_ptr<LinSolver>(
                    new KSPSolver(solverName, configFile));
            break;
        case petibm::utilities::types::ExecuteType::GPU:
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

} // end of namespace linsolvers
} // end of namespace petibm
