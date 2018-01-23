/** 
 * \file linsolver.cpp
 * \brief Implementations of LinSolverBase, LinSolver, and factory function.
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \copyright MIT.
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

// implement LinSolverBase::destroy
PetscErrorCode LinSolverBase::destroy()
{
    PetscFunctionBeginUser;
    name = config = type = "";
    PetscFunctionReturn(0);
} // destroy
    
    
// implement LinSolverBase::printInfo
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
} // printInfo


// implement LinSolverBase::getType
PetscErrorCode LinSolverBase::getType(std::string &_type) const
{
    PetscFunctionBeginUser;
    _type = type;
    PetscFunctionReturn(0);
} // getType


// implement petibm::linsolver::createLinSolver
PetscErrorCode createLinSolver(const std::string &solverName, 
        const YAML::Node &node, type::LinSolver &solver)
{
    PetscFunctionBeginUser;
    
    std::string     key, config, type;
    
    key = solverName + "Solver";
    
    // set up type
    type = node["parameters"][key]["type"].as<std::string>("CPU");
        
    // set up the path to config file
    config = node["parameters"][key]["config"].as<std::string>("None");
    if (config[0] != '/' && config != "None")
    	config = node["directory"].as<std::string>() + "/" + config;

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
} // createLinSolver

} // end of namespace linsolver
} // end of namespace petibm
