/**
 * \file solution.h
 * \brief Definition of the class petibm::solution::SolutionBase,
 *        the type definition petibm::type::Solution,
 *        and the factory function petibm::solution::createSolution.
 * function. \copyright Copyright (c) 2016-2018, Barba group. All rights
 * reserved. \license BSD 3-Clause License.
 */

#pragma once

#include <memory>
#include <string>

#include <petscsys.h>

#include <petibm/mesh.h>

/**
 * \defgroup solutionModule Solution holders
 * \brief Simple structures to hold flow field solutions.
 *
 * Objects hold the solution of the pressure scalar field
 * and the velocity vector field.
 * These flow fields are stored in distributed PETSc Vec objects.
 * Objects also handles I/O of the flow field solutions.
 *
 * petibm::solution::SolutionSimple is currently the only adaptation
 * of the abstract class petibm::solution::SolutionBase.
 * API users should use the type definition petibm::solution::Solution
 * to create a solution object through the factory function
 * petibm::solution::createSolution.
 * Flow solvers implemented in PetIBM currently do not need other
 * functionalities than those available in the abstract class and
 * its only adaptation.
 * We use an abstract class for the potential of expanding this category
 * in the future.
 *
 * \see petibm::type::Solution, petibm::solution::createSolution
 * \ingroup petibm
 */

namespace petibm
{
/**
 * \brief Collection of classes and functions regarding solution holders.
 *
 * \see solutionModule, petibm::type::Solution, petibm::solution::createSolution
 * \ingroup solutionModule
 */
namespace solution
{
/**
 * \brief Base (abstract) class for different solution holders.
 *
 * \see solutionModule, petibm::type::Solution, petibm::solution::createSolution
 * \ingroup solutionModule
 */
class SolutionBase
{
public:
    /** \brief Number of dimensions. */
    PetscInt dim;

    /** \brief Packed PETSc Vec object for the velocity vector field. */
    Vec UGlobal;

    /** \brief PETSc Vec for the pressure scalar field. */
    Vec pGlobal;

    /** \brief String containing information about the solution. */
    std::string info;

    /** \brief Default constructor. */
    SolutionBase() = default;

    /**
     * \brief Constructor using a Cartesian mesh.
     *
     * \param mesh [in] Cartesian mesh object.
     */
    SolutionBase(const type::Mesh &mesh){};

    /** \brief Default destructor. */
    virtual ~SolutionBase();

    /**
     * \brief Manually destroy data.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode destroy();

    /**
     * \brief Print information about the solution to standard output.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode printInfo() const;

    /**
     * \brief Set initial conditions of the flow fields.
     *
     * \param node [in] YAML node with flow settings.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode setInitialConditions(const YAML::Node &node) = 0;

    /**
     * \brief Convert velocity fluxes to velocity components.
     * 
     * It assumes that UGlobal currently contains velocity fluxes
     * and converts (inplace) to velocity components.
     *
     * \param Rinv [in] Operator to convert fluxes to components.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode convert2Velocity(const Mat &Rinv) = 0;

    /**
     * \brief Convert velocity components to velocity fluxes.
     * 
     * It assumes that UGlobal currently contains velocity components
     * and converts (inplace) to velocity fluxes.
     *
     * \param R [in] Operator to convert components to fluxes.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode convert2Flux(const Mat &R) = 0;

    /**
     * \brief Write flow field solutions to a file.
     *
     * Currently only supports HDF5 format.
     *
     * \param filePath [in] Path of the file to write in.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode write(const std::string &filePath) const = 0;

    /**
     * \brief Read the flow field solutions from a file.
     *
     * Currently only supports HDF5 format.
     *
     * \param filePath [in] Path of the file to read from.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode read(const std::string &filePath) = 0;

protected:
    /**
     * \brief Initialize the flow field solutions.
     *
     * \param mesh [in] Cartesian mesh object.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode init(const type::Mesh &mesh) = 0;

    /** \brief MPI communicator. */
    MPI_Comm comm;

    /** \brief Size of MPI communicator. */
    PetscMPIInt mpiSize;

    /** \brief Rank of the local process. */
    PetscMPIInt mpiRank;

    /** \brief Shared pointer to the underlying Cartesian mesh object. */
    type::Mesh mesh;

};  // SolutionBase
}  // end of namespace solution

namespace type
{
/**
 * \brief Type definition of solution object.
 *
 * Users should use the function petibm::solution::createSolution
 * to create a petibm::solution::Solution object.
 *
 * Example usage:
 * \code
 * PetscErrorCode ierr;
 * YAML::Node config;
 * petibm::type::Mesh mesh;
 * petibm::type::Solution solution;
 *
 * // get the simulation settings into YAML node config
 * // create Mesh with petibm::mesh::createMesh
 *
 * // create a Solution instance
 * ierr = petibm::solution::createSolution(mesh, solution); CHKERRQ(ierr);
 * // set initial conditions
 * ierr = solution->setInitialConditions(config); CHKERRQ(ierr);
 * // write the flow field solutions into a file test.h5
 * ierr = solution->write("./test.h5"); CHKERRQ(ierr);
 * \endcode
 *
 * \see solutionModule, petibm::solution::SolutionBase,
 * petibm::solution::createSolution \ingroup solutionModule
 */
typedef std::shared_ptr<solution::SolutionBase> Solution;
}  // end of namespace type

namespace solution
{
/**
 * \brief Factory function to create a petibm::solution::Solution object.
 *
 * The function create a solution object based on the provided underlying
 * Cartesian mesh.
 * 
 * \param mesh [in] Cartesian mesh object.
 * \param solution [out] Solution object.
 *
 * \return PetscErrorCode.
 *
 * \see solutionModule, petibm::type::Solution
 * \ingroup solutionModule
 */
PetscErrorCode createSolution(const type::Mesh &mesh, type::Solution &solution);
}  // end of namespace solution

}  // end of namespace petibm
