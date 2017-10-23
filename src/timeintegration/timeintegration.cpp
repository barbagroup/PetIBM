/***************************************************************************//**
 * \file timeintegration.h
 * \author Anush Krishnan (anush@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Definition TimeIntegration related code.
 */


// here goes headers from our PetIBM
# include <petibm/timeintegration.h>


namespace petibm
{
namespace timeintegration
{
PetscErrorCode TimeIntegrationBase::printInfo() const
{
    PetscFunctionBeginUser;
    
    PetscErrorCode ierr;
    
    std::string info;
    
    info += (std::string(80, '=') + "\n");
    info += ("Time Integration [" + name + "]\n");
    info += (std::string(80, '=') + "\n");
    
    info += ("\tScheme: " + scheme + "\n\n");
    info += ("\tCoefficient of implicit term: " + 
            std::to_string(implicitCoeff) + "\n\n");
    info += "\tCoefficients of Explicit terms: [";
    for(auto it: explicitCoeffs) info += (std::to_string(it) + ", ");
    info += "]\n\n";
    
    ierr = PetscPrintf(PETSC_COMM_WORLD, "%s", info.c_str()); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

PetscErrorCode createTimeIntegration(const std::string &name,
        const YAML::Node &node, type::TimeIntegration &integration)
{
    PetscFunctionBeginUser;

    std::string     scheme;
    
    if (! node["parameters"].IsDefined())
    {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                "Could not find the key \"parameters\" in the YAML node "
                "passed to the function \"createTimeIntegration\"!");
    }
    
    if (! node["parameters"][name].IsDefined())
    {
        SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                "Could not find the key \"%s\" under sub-node \"parameters\" "
                "in the YAML node passed to the function "
                "\"createTimeIntegration\"!", name.c_str());
    }
    
    scheme = node["parameters"][name].as<std::string>();

    // set up coefficients.
    if (scheme == "EULER_EXPLICIT")
        integration = std::make_shared<Euler_Explicit>(name);
    else if (scheme == "EULER_IMPLICIT")
        integration = std::make_shared<Euler_Implicit>(name);
    else if (scheme == "ADAMS_BASHFORTH_2")
        integration = std::make_shared<Adams_Bashforth_2>(name);
    else if (scheme == "CRANK_NICOLSON")
        integration = std::make_shared<Crank_Nicolson>(name);
    else
    {
        SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                "The time integration scheme \"%s\" does not exist.\n",
                scheme.c_str());
    }


    PetscFunctionReturn(0);
}

} // end of namespace timeintegration
} // end of namespace petibm
