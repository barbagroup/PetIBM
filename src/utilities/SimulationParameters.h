/***************************************************************************//**
 * \file SimulationParameters.h
 * \author Anush Krishnan (anush@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Definition of the class `SimulationParameters`.
 */


# pragma once

// here goes C++ STL
# include <iostream>
# include <string>
# include <vector>
# include <experimental/filesystem>

// here goes PETSc headers
# include <petscsys.h>

// here goes YAML header
# include <yaml-cpp/yaml.h>

// here goes headers from our PetIBM
# include "types.h"


/**
 * \class SimulationParameters
 * \brief Stores various parameters used in the simulation.
 */
class SimulationParameters
{
public:

    /** \brief directory of the simulation. */
    stdfs::path                 caseDir;


    /** \brief output settings. */
    types::OutputInfo       output;


    /** \brief info of velocity solver. */
    types::LinSolverInfo    vSolver;

    /** \brief info of poisson solver. */
    types::LinSolverInfo    pSolver;



    /** \brief info of numerical schemes. */
    types::SchemeInfo       schemes;


    /** \brief time increment. */
    types::SteppingInfo     step;


    /** \brief a string of information for this object. */
    std::string             info;
  

    /** \brief Default constructor. */
    SimulationParameters();

    /**
     * \brief constructor with given config YAML file.
     *
     * \param world MPI communicator.
     * \param config YAML config file.
     */
    SimulationParameters(const MPI_Comm &world, const std::string &config);

    /**
     * \brief constructor with given YAML node.
     *
     * \param world MPI communicator.
     * \param config YAML config node.
     * \param dir the path to the case directory; if empty, then there must be 
     *            a key "caseDir" in the YAML Node and at the same level as the
     *            key "simulationParameters".
     */
    SimulationParameters(const MPI_Comm &world, 
            const YAML::Node &config, const std::string &dir="");

    /** \brief destructor */
    ~SimulationParameters();


    /**
     * \brief initialization of the SimulationParameter object
     *
     * \param world MPI communicator.
     * \param config YAML config node.
     * \param dir the path to the case directory.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode init(const MPI_Comm &world, 
            const YAML::Node &config, const std::string &dir="");


    /**
     * \brief print information using only the master process.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode printInfo() const;

protected:

    std::shared_ptr<const MPI_Comm>         comm;

    PetscMPIInt                             mpiSize,
                                            mpiRank;

    /**
     * \brief create a std::string holding the information for output.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode createInfoString();

    /**
     * \brief check if the setting use HDF5 and if we were compiled with HDF5.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode checkHDF5();

    /**
     * \brief check if the setting use GPU and if we were compiled with AmgX.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode checkAmgX();

}; // SimulationParameters


/**
 * \brief overloading of output stream operator `<<`.
 *
 * \param os the output stream.
 * \param param a SimulationParameters onject.
 *
 * \return 
 */
std::ostream &operator<<(std::ostream &os, const SimulationParameters &param);
