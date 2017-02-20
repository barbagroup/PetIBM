/***************************************************************************//**
 * \file FlowDescription.cpp
 * \author Anush Krishnan (anush@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Implementation of the methods of the class `FlowDescription`.
 */

// here goes C++ STL
# include <sstream>

// here goes our own headers
# include "FlowDescription.h"
# include "parser.h"


/** \copydoc FlowDescription::FlowDescription() */
FlowDescription::FlowDescription() { }


/** \copydoc FlowDescription::FlowDescription(const std::string &) */
FlowDescription::FlowDescription(const std::string &YAMLfile)
{
    YAML::Node      node = YAML::LoadFile(YAMLfile);

    // see if the key "flowDescription" is used in the YAML file
    if (node["flowDescription"].IsDefined())
        // not sure why GCC can't understand node["flowDescription"] is a node
        FlowDescription(YAML::Node(node["flowDescription"]));
    else // pass the whole node
        FlowDescription(node);
}


/** \copydoc FlowDescription::FlowDescription(const YAML::Node &) */
FlowDescription::FlowDescription(const YAML::Node &node)
{
    PetscPrintf(PETSC_COMM_WORLD, "Creating a FlowDescription object ... ");

    if (node["flowDescription"].IsDefined())
        parser::parseFlowDescription(node["flowDescription"], 
                dim, nu, customIC, IC, pertb, nBCs, BCInfo);
    else
        parser::parseFlowDescription(node, 
                dim, nu, customIC, IC, pertb, nBCs, BCInfo);

    checkPeriodicity();
    createInfoString();

    PetscPrintf(PETSC_COMM_WORLD, "done.\n");
}


/** \copydoc FlowDescription::~FlowDescription() */
FlowDescription::~FlowDescription() {}


/** \copydoc FlowDescription::checkPeriodicity() */
PetscErrorCode FlowDescription::checkPeriodicity()
{
    using namespace types;

    PetscFunctionBeginUser;

    // if a component at one boundary is periodic, then other components should
    // also be periodic
    for(auto loc: BCInfo)
    {
        int     pCount = 0;

        // count the number of periodic BCs at current boundary
        for(auto comp: loc.second)
            if (comp.second.type == BCType::PERIODIC) pCount += 1;

        // check if all components are all periodic or all not
        if ((pCount != 0) && (pCount != dim))
            SETERRQ1(PETSC_COMM_WORLD, 60, "Not all BC types on boundary %s "
                    "are periodic !", bl2str[loc.first].c_str());
    }

    // if one boundary is periodic, then the counterpart should also be
    BCLoc   i, j;
    for(i=BCLoc::XMINUS, j=BCLoc::XPLUS; i<nBCs; i=BCLoc(int(i)+2), j=BCLoc(int(j)+2))
    {
        if ((BCInfo[i][u].type == BCType::PERIODIC) &&
                (BCInfo[j][u].type != types::BCType::PERIODIC))
            SETERRQ2(PETSC_COMM_WORLD, 60, "%s is a periodic boundary, while "
                    "its counterpart %s is not !!", 
                    types::bl2str[i].c_str(), types::bl2str[j].c_str());

        if ((BCInfo[j][u].type == BCType::PERIODIC) &&
                (BCInfo[i][u].type != types::BCType::PERIODIC))
            SETERRQ2(PETSC_COMM_WORLD, 60, "%s is a periodic boundary, while "
                    "its counterpart %s is not !!", 
                    types::bl2str[j].c_str(), types::bl2str[i].c_str());
    }

    PetscFunctionReturn(0);
}


/** \copydoc FlowDescription::createInfoString() */
PetscErrorCode FlowDescription::createInfoString()
{
    using namespace types;

    PetscFunctionBeginUser;

    std::stringstream       ss;

    ss << "Flow Descriptions: " << std::endl;
    ss << "\tDimension: " << dim << std::endl;
    ss << "\tViscosity: " << nu << std::endl;
    ss << "\tCustomized IC: " << (customIC? "true" : "false") << std::endl;
    ss << "\tIC: [ ";
    for(auto it: IC) ss << it << " ";
    ss << "]" << std::endl;
    ss << "\tPerturbation:" << std::endl;
    ss << "\t\tFrequency: " << pertb.freq << std::endl;
    ss << "\t\tAmplitude: " << pertb.amp << std::endl;
    ss << "\tBoundary conditions:" << std::endl;
    for(auto loc: BCInfo)
    {
        ss << "\t\t" << bl2str[loc.first] << ":" << std::endl;
        for(auto comp: loc.second)
            ss << "\t\t\t" << fd2str[comp.first] << ": [ "
                << bt2str[comp.second.type] << ", "
                << comp.second.value << " ]" << std::endl;
    }

    info = ss.str();

    PetscFunctionReturn(0);
}


/** \copydoc FlowDescription::printInfo() */
PetscErrorCode FlowDescription::printInfo()
{
    PetscFunctionBeginUser;

    PetscErrorCode  ierr;

    ierr = PetscPrintf(PETSC_COMM_WORLD, info.c_str()); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/** \copydoc operator<<(std::ostream &, const FlowDescription &) */
std::ostream &operator<<(std::ostream &os, const FlowDescription &flow)
{
    os << flow.info;
    return os;
}
