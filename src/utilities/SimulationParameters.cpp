/***************************************************************************//**
 * \file SimulationParameters.cpp
 * \author Anush Krishnan (anush@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Implementation of the methods of the class `SimulationParameters`.
 */


// here goes headers from our PetIBM
# include "SimulationParameters.h"
# include "parser.h"


using namespace types;


/** \copydoc SimulationParameters::SimulationParameters() */
SimulationParameters::SimulationParameters() = default;


/** \copydoc SimulationParameters::SimulationParameters(
 * const MPI_Comm &world, const YAML::Node &config, const std::string &dir). */
SimulationParameters::SimulationParameters(const MPI_Comm &world, 
        const YAML::Node &config, const std::string &dir)
{
    init(world, config, dir);
}


/** \copydoc SimulationParameters::~SimulationParameters() */
SimulationParameters::~SimulationParameters() = default;


/** \copydoc SimulationParameters::init() */
PetscErrorCode SimulationParameters::init(const MPI_Comm &world, 
        const YAML::Node &node, const std::string &dir)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    // store the address of the communicator
    // note: this is a bad practice; shared_ptr is not for stack variables!!
    comm = std::shared_ptr<const MPI_Comm>(&world, [](const MPI_Comm*){});

    // set rank and size
    ierr = MPI_Comm_size(*comm, &mpiSize); CHKERRQ(ierr);
    ierr = MPI_Comm_rank(*comm, &mpiRank); CHKERRQ(ierr);
    
    // store the case directory
    caseDir = stdfs::path(dir);

    // get data from the YAML node
    ierr = parser::parseSimulationParameters(
            node, output, vSolver, pSolver, schemes, step); CHKERRQ(ierr);

    // check HDF5
    ierr = checkHDF5(); CHKERRQ(ierr);

    // check AmgX
    ierr = checkAmgX(); CHKERRQ(ierr);

    // create a string for printing or output stream
    ierr = createInfoString(); CHKERRQ(ierr);

    ierr = MPI_Barrier(*comm); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/** \copydoc SimulationParameters::checkHDF5(). */
PetscErrorCode SimulationParameters::checkHDF5()
{
    PetscFunctionBeginUser;

#ifndef PETSC_HAVE_HDF5
    if (output.format == OutputFormat::HDF5)
        SETERRQ(PETSC_COMM_WORLD, 56, "PETSc was not built with HDF5; "
                "you cannot use `outputFormat: hdf5`");
#endif

    PetscFunctionReturn(0);
}


/** \copydoc SimulationParameters::checkAmgX(). */
PetscErrorCode SimulationParameters::checkAmgX()
{
    PetscFunctionBeginUser;

#ifndef HAVE_AMGX
    if ((pSolver.type == ExecuteType::GPU) ||
            (vSolver.type == ExecuteType::GPU))
        SETERRQ(PETSC_COMM_WORLD, 56, "PetIBM was not built with AmgX enabled; "
                "you cannot use GPU linear solvers");
#endif

    PetscFunctionReturn(0);
}


/** \copydoc SimulationParameters::createInfoString(). */
PetscErrorCode SimulationParameters::createInfoString()
{
    PetscFunctionBeginUser;

    std::stringstream       ss;


    ss << std::string(80, '=') << std::endl;
    ss << "Simulation Parameters:" << std::endl;
    ss << std::string(80, '=') << std::endl;

    ss << "\tCase Directory:" << std::endl;
    ss << "\t\t" << caseDir << std::endl;
    ss << std::endl;

    ss << "\tFlow solver: " << ibm2str[schemes.ibm] << std::endl;
    ss << std::endl;

    ss << "\tTemporal Discretization: " << std::endl;
    ss << "\t\tConvection: " << ts2str[schemes.convection] << std::endl;
    ss << "\t\tDiffusion: " << ts2str[schemes.diffusion] << std::endl;
    ss << std::endl;

    ss << "\tTime-Stepping Control: " << std::endl;
    ss << "\t\tTime Step: " << step.dt << std::endl;
    ss << "\t\tStarting Step: " << step.nStart << std::endl;
    ss << "\t\tTotal Time Steps: " << step.nTotal << std::endl;
    ss << "\t\tNumber of Steps to Save Solutions: " << step.nSave << std::endl;
    ss << "\t\tNumber of Steps to Save Restart Data: " << step.nRestart << std::endl;
    ss << std::endl;

    ss << "\tVelocity Linear Solver:" << std::endl;
    ss << "\t\tType: " << et2str[vSolver.type] << std::endl;
    ss << "\t\tConfiguration File: " << vSolver.config << std::endl;
    ss << std::endl;

    ss << "\tPoisson Linear Solver:" << std::endl;
    ss << "\t\tType: " << et2str[pSolver.type] << std::endl;
    ss << "\t\tConfiguration File: " << pSolver.config << std::endl;
    ss << std::endl;

    ss << "\tOutput Format Control:" << std::endl;
    ss << "\t\tFormat: " << out2str[output.format] << std::endl;
    ss << "\t\tOutput Flux: " 
        << (output.outputFlux ? "true" : "false") << std::endl;
    ss << "\t\tOutput Velocity: " 
        << (output.outputVelocity ? "true" : "false") << std::endl;
    ss << std::endl;

    info = ss.str();

    PetscFunctionReturn(0);
}


/** \copydoc SimulationParameters::printInfo(). */
PetscErrorCode SimulationParameters::printInfo() const
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    ierr = PetscPrintf(*comm, info.c_str()); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/** \copyfoc std::ostream &operator<<(std::ostream &, const SimulationParameters &). */
std::ostream &operator<<(std::ostream &os, const SimulationParameters &param)
{
    os << param.info;
    return os;
}
