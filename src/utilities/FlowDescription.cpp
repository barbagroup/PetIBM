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


using namespace types;


/** \copydoc FlowDescription::FlowDescription() */
FlowDescription::FlowDescription() { }


/** \copydoc FlowDescription::FlowDescription(const MPI_Comm &, const std::string &) */
FlowDescription::FlowDescription(const MPI_Comm &world, const std::string &YAMLfile)
{
    YAML::Node      node = YAML::LoadFile(YAMLfile);

    FlowDescription(world, node);
}


/** \copydoc FlowDescription::FlowDescription(const MPI_Comm &, const YAML::Node &) */
FlowDescription::FlowDescription(const MPI_Comm &world, const YAML::Node &node)
{
    if (node["flowDescription"].IsDefined())
        init(world, node["flowDescription"]);
    else
        init(world, node);
}


/** \copydoc FlowDescription::init(const MPI_Comm &, const YAML::Node &) */
PetscErrorCode FlowDescription::init(const MPI_Comm &world, const YAML::Node &node)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    // store the address of the communicator
    // note: this is a bad practice; shared_ptr is not for stack variables!!
    comm = std::shared_ptr<const MPI_Comm>(&world, [](const MPI_Comm*){});

    // set rank and size
    ierr = MPI_Comm_size(*comm, &mpiSize); CHKERRQ(ierr);
    ierr = MPI_Comm_rank(*comm, &mpiRank); CHKERRQ(ierr);

    // parse YAML node
    parser::parseFlowDescription(node, 
            dim, nu, customIC, IC, pertb, nBCs, BCInfo);

    // check if the setting of periodic BCs is correct
    ierr = checkPeriodicity(); CHKERRQ(ierr);

    // create a string to hold information
    ierr = createInfoString(); CHKERRQ(ierr);

    ierr = MPI_Barrier(*comm); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/** \copydoc FlowDescription::~FlowDescription() */
FlowDescription::~FlowDescription() {}


/** \copydoc FlowDescription::checkPeriodicity() */
PetscErrorCode FlowDescription::checkPeriodicity()
{
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
                (BCInfo[j][u].type != BCType::PERIODIC))
            SETERRQ2(PETSC_COMM_WORLD, 60, "%s is a periodic boundary, while "
                    "its counterpart %s is not !!", 
                    bl2str[i].c_str(), bl2str[j].c_str());

        if ((BCInfo[j][u].type == BCType::PERIODIC) &&
                (BCInfo[i][u].type != BCType::PERIODIC))
            SETERRQ2(PETSC_COMM_WORLD, 60, "%s is a periodic boundary, while "
                    "its counterpart %s is not !!", 
                    bl2str[j].c_str(), bl2str[i].c_str());
    }

    PetscFunctionReturn(0);
}


/** \copydoc FlowDescription::createInfoString() */
PetscErrorCode FlowDescription::createInfoString()
{
    PetscFunctionBeginUser;

    std::stringstream       ss;


    ss << std::string(80, '=') << std::endl;
    ss << "Flow Descriptions: " << std::endl;
    ss << std::string(80, '=') << std::endl;

    ss << "\tDimension: " << dim << std::endl;
    ss << std::endl;

    ss << "\tViscosity: " << nu << std::endl;
    ss << std::endl;

    ss << "\tCustomized IC: " << (customIC? "true" : "false") << std::endl;
    ss << "\tIC: [ ";
    for(auto it: IC) ss << it << " ";
    ss << "]" << std::endl;
    ss << std::endl;

    ss << "\tPerturbation:" << std::endl;
    ss << "\t\tFrequency: " << pertb.freq << std::endl;
    ss << "\t\tAmplitude: " << pertb.amp << std::endl;
    ss << std::endl;

    ss << "\tBoundary conditions:" << std::endl;
    for(auto loc: BCInfo)
    {
        ss << "\t\t" << bl2str[loc.first] << ":" << std::endl;
        for(auto comp: loc.second)
            ss << "\t\t\t" << fd2str[comp.first] << ": [ "
                << bt2str[comp.second.type] << ", "
                << comp.second.value << " ]" << std::endl;
    }
    ss << std::endl;

    info = ss.str();

    PetscFunctionReturn(0);
}


/** \copydoc FlowDescription::printInfo() */
PetscErrorCode FlowDescription::printInfo() const
{
    PetscFunctionBeginUser;

    PetscErrorCode  ierr;

    ierr = PetscPrintf(*comm, info.c_str()); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/** \copydoc operator<<(std::ostream &, const FlowDescription &) */
std::ostream &operator<<(std::ostream &os, const FlowDescription &flow)
{
    os << flow.info;
    return os;
}
