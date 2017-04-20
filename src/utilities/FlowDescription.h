/***************************************************************************//**
 * \file FlowDescription.h
 * \author Anush Krishnan (anush@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Definition of the class `FlowDescription`.
 */

# pragma once

// here goes C++ STL
# include <iostream>
# include <string>

// here goes PETSc
# include <petscsys.h>

// here goes yaml-cpp
# include <yaml-cpp/yaml.h>

// here goes our headers
# include "types.h"


/**
 * \class FlowDescription
 * \brief Stores information that describes the flow.
 *
 * \todo remember to check if the capability of specifying a different 
 * flowDescription.yaml is implemented in other place. 
 *
 * \todo think about whether we should use PetscInt to store BCinfo or use
 * types::BCLoc and types::Field. Using PetscInt allow us to use
 * std::vector, whil using types::BCLoc and types::VeclocityComponent will 
 * require std::map.
 */
class FlowDescription
{
public:

    /** \brief dimension of the flow. */
    PetscInt                    dim;

    /** \brief kinematic viscosity. */
    PetscReal                   nu;

    /** \brief flag indicating if IC will be customized (not supported yet). */
    PetscBool                   customIC;

    /** \brief initial conditions. */
    types::RealVec1D            IC;

    /** \brief perturbation, including frequency and amplitude. */
    types::Perturbation         pertb;

    /** \brief number of boundaries, should be dim * 2. */
    PetscInt                    nBCs;

    /** \brief information of BCs obtained from YAML file. */
    types:: BCInfoHolder        BCInfo;

    /** \brief a string or flow information that can be used for printing. */
    std::string                 info;


    /** \brief default constructor. */
    FlowDescription();


    /** \brief constructor for using YAML node.
     *
     * \ details Either a flowDescription node or a YAML node that has a key 
     * `flowDescription` is accepted.
     *
     * \param world MPI communicator.
     * \param node a YAML node.
     */
    FlowDescription(const MPI_Comm &world, const YAML::Node &node);

    /** \brief destructor. */
    ~FlowDescription();

    /** \brief initialization.
     *
     * \param world MPI communicator.
     * \param node a YAML node
     */
    PetscErrorCode init(const MPI_Comm &world, const YAML::Node &node);

    /** \brief print information using PETSc printf.
     *
     * \details This is for parallel simulations. With PETSc printf, only the 
     * master node will print informations.
     *
     * \return PetscErrorCode
     */
    PetscErrorCode printInfo() const;

protected:

    std::shared_ptr<const MPI_Comm>         comm;

    PetscMPIInt                             mpiSize,
                                            mpiRank;

    /** \brief check if users' settings of periodicity are correct.
     *
     * \return PetscErrorCode
     */
    PetscErrorCode checkPeriodicity();

    /** \brief create a std::string for information.
     *
     * \details This string then can be used in both PETSc printf or C++ 
     * output streams.
     *
     * \return PetscErrorCode
     */
    PetscErrorCode createInfoString();

}; // FlowDescription


/** \brief output stream support for FlowDescription onjects.
 *
 * \details The content sent to output stream will be the info string.
 *
 * \param os stdout or other output streams.
 * \param flow a FlowDescription instance.
 *
 * \return output stream
 */
std::ostream &operator<<(std::ostream &os, const FlowDescription &flow);
