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
    
    
PetscErrorCode LinSolverBase::printInfo() const
{
    PetscFunctionBeginUser;
    
    PetscErrorCode  ierr;
    
    std::string info = "";
    info += (std::string(80, '=') + "\n");
    info += ("Linear Solver " + name + ":\n");
    info += (std::string(80, '=') + "\n");
    
    info += ("\tType: " + type + "\n\n");
    info += ("\tConfig file: " + config + "\n\n");
    
    ierr = PetscPrintf(PETSC_COMM_WORLD, "%s", info.c_str()); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}


// factory function for creating a LinSolver
PetscErrorCode createLinSolver(const std::string &solverName, 
        const YAML::Node &node, type::LinSolver &solver)
{
    PetscFunctionBeginUser;
    
    std::string     key, config, type;
    
    key = solverName + "Solver";
    
    if (! node["parameters"].IsDefined())
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                "Could not find the key \"parameters\" in the YAML node passed "
                "to the function \"createLinSolver\"\n");
    
    if (! node["parameters"][key].IsDefined())
        SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                "Could not find the key \"%s\" under \"parameters\" in the YAML "
                "node passed to the function \"createLinSolver\"\n", key.c_str());
    
    if (! node["parameters"][key]["type"].IsDefined())
        SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                "Could not find the key \"type\" under the settings for linear "
                "solver \"%s\"\n", solverName.c_str());
    
    if (! node["parameters"][key]["config"].IsDefined())
        SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                "Could not find the key \"config\" under the settings for linear "
                "solver \"%s\"\n", solverName.c_str());
    
    // set up type
    type = node["parameters"][key]["type"].as<std::string>();
        
    // set up the path to config file
    config = node["parameters"][key]["config"].as<std::string>();
    if (config[0] != '/') config = node["directory"].as<std::string>() + "/" + config;

    // factory
    if (type == "CPU")
        solver = std::make_shared<LinSolverKSP>(solverName, config);
    else if (type == "GPU")
# ifdef HAVE_AMGX
        solver = std::make_shared<LinSolverAmgX>(solverName, config);
# else
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                "AmgX solver is used, while PetIBM is not compiled with AmgX.");
# endif
    else
        SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                "Unrecognized value \"%s\" of the type of the linear solver "
                "\"%s\"\n", type.c_str(), solverName.c_str());

    PetscFunctionReturn(0);
}

} // end of namespace linsolver
} // end of namespace petibm