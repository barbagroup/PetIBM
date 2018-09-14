/**
 * \file solution.h
 * \brief Definition of solution::SolutionBase, type::Solution, and factory function.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 */


# pragma once

// C++ STL
# include <string>
# include <memory>

// PETSc
# include <petscsys.h>

// PetIBM
# include <petibm/mesh.h>


/**
 * \defgroup solutionModule Solution holders
 * \brief Simple and useful structures for holding solutions.
 *
 * Useful objects that hold solution vectors. These objects also handle the I/O
 * of solutions.
 * 
 * Though there is only one implementation in this category -- 
 * petibm::solution::SolutionSimple, API users should still use 
 * petibm::type::Solution to access this only-one implementation.
 * The flow solvers implemented in PetIBM
 * currently don't need other types of implementation, but using an abstract 
 * class provides potential to expand this category in the future.
 * 
 * \see petibm::type::Solution, petibm::solution::createSolution
 * \ingroup petibm
 */


namespace petibm
{
/** 
 * \brief Collection of classes and functions regarding solution holders.
 * \see solutionModule, petibm::type::Solution, petibm::solution::createSolution
 * \ingroup solutionModule
 */
namespace solution
{
/**
 * \brief Base (abstract) class for different solution holders.
 * \see solutionModule, petibm::type::Solution, petibm::solution::createSolution
 * \ingroup solutionModule
 */
class SolutionBase
{
public:

    /** \brief Dimension. */
    PetscInt                dim;

    /** \brief Packed PETSc Vec for velocities. */
    Vec                     UGlobal;

    /** \brief A PETSc Vec for pressure field. */
    Vec                     pGlobal;

    /** \brief A std::string containing information. */
    std::string             info;
    
    
    /** \brief Default constructor. */
    SolutionBase() = default;

    /**
     * \brief Construct through a type::Mesh.
     *
     * \param mesh [in] a type::Mesh object.
     */
    SolutionBase(const type::Mesh &mesh) {};

    /** \brief Default destructor. */
    virtual ~SolutionBase();

    /**
     * \brief Manually destroy data.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode destroy();
    
    /**
     * \brief Print information to standard output.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode printInfo() const;

    /**
     * \brief Apply ICs through the settings in a YAML node.
     *
     * \param node [in] a YAML::Node containing settings of the flow.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode applyIC(const YAML::Node &node) = 0;

    /**
     * \brief Assume the values in UGlobal are fluxes and convert them to 
     *        velocities.
     *
     * \param Rinv [in] a PETSc Mat representing converting operator.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode convert2Velocity(const Mat &Rinv) = 0;

    /**
     * \brief Assume the values in UGlobal are velocities and convert them to 
     *        fluxes.
     *
     * \param R [in] a PETSc Mat representing converting operator.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode convert2Flux(const Mat &R) = 0;
    
    /**
     * \brief Write solution vectors to a given file.
     * 
     * Now only support HDF5 format.
     *
     * \param filePath [in] path of the file to write in.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode write(const std::string &filePath) const = 0;
    
    /**
     * \brief Read solution vectors from a given file.
     * 
     * Now only support HDF5 format.
     *
     * \param filePath [in] path of the file to read from.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode read(const std::string &filePath) = 0;


protected:

    /**
     * \brief Underlying initialization function.
     *
     * \param mesh [in] a type::Mesh.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode init(const type::Mesh &mesh) = 0;

    /** \brief MPI communicator. */
    MPI_Comm                                comm;

    /** \brief Size of MPI communicator. */
    PetscMPIInt                             mpiSize;
    
    /** \brief The rank of this process. */
    PetscMPIInt                             mpiRank;

    /** \brief A std::shared_ptr to underlying mesh. */
    type::Mesh                              mesh;

}; // SolutionBase
} // end of namespace solution


namespace type
{
    /**
     * \brief Type definition of Solution.
     * 
     * Please use petibm::solution::createSolution to create a Solution object.
     * 
     * Example usage:
     * \code
     * PetscErrorCode ierr;
     * petibm::type::Mesh mesh;
     * petibm::type::Solution soln;
     * 
     * // create Mesh with petibm::mesh::createMesh
     * 
     * // create a Solution instance
     * ierr = petibm::solution::createSolution(mesh, soln); CHKERRQ(ierr);
     * ierr = soln->applyIC(config); CHKERRQ(ierr); // config is a YAML::Node
     * ierr = soln->write("./test"); CHKERRQ(ierr); // write a HDF5 file, ./test.h5
     * 
     * \endcode
     * 
     * \see solutionModule, petibm::solution::SolutionBase, petibm::solution::createSolution
     * \ingroup solutionModule
     */
    typedef std::shared_ptr<solution::SolutionBase> Solution;
} // end of namespace type


namespace solution
{
    /**
     * \brief A factory function for creating Solution objects.
     * \param mesh [in] a type::Mesh.
     * \param solution [out] resulting Solution object.
     * \return PetscErrorCode.
     * 
     * This function will create a Solution instance based on the information in
     * the provided mesh object.
     * 
     * \see solutionModule, petibm::type::Solution
     * \ingroup solutionModule
     */
    PetscErrorCode createSolution(
            const type::Mesh &mesh, type::Solution &solution);
} // end of namespace solution

} // end of namespace petibm
