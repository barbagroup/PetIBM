/***************************************************************************//**
 * \file solutions.h
 * \author Anush Krishnan (anush@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Definition of the class `Solutions`.
 */


# pragma once

// C++ STL
# include <string>
# include <memory>

// PETSc
# include <petscsys.h>

// PetIBM
# include <petibm/mesh.h>


namespace petibm
{
namespace solution
{
/** \brief base abstract class for solutions. */
class SolutionBase
{
public:

    /** \brief dimension. */
    PetscInt                dim;

    /** \brief packed PETSc Vec for velocities. */
    Vec                     UGlobal;

    /** \brief a PETSc Vec for pressure field. */
    Vec                     pGlobal;

    /** \brief a std::string containing information. */
    std::string             info;
    
    
    /** \brief default constructor. */
    SolutionBase() = default;

    /**
     * \brief construct through a type::Mesh.
     *
     * \param mesh [in] a type::Mesh object.
     */
    SolutionBase(const type::Mesh &mesh) {};

    /** \brief default destructor. */
    virtual ~SolutionBase() = default;
    
    /**
     * \brief print information to standard output.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode printInfo() const;

    /**
     * \brief apply ICs through the settings in a YAML node.
     *
     * \param node [in] a YAML::Node.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode applyIC(const YAML::Node &node) = 0;

    /**
     * \brief assume the values in UGlobal are fluxes and convert them to 
     *        velocities.
     *
     * \param Rinv [in] a PETSc Mat representing converting operator.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode convert2Velocity(const Mat &Rinv) = 0;

    /**
     * \brief assume the values in UGlobal are velocities and convert them to 
     *        fluxes.
     *
     * \param R [in] a PETSc Mat representing converting operator.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode convert2Flux(const Mat &R) = 0;
    
    /**
     * \brief write solution vectors to a given file.
     * 
     * Now only support HDF5 format.
     *
     * \param file [in] the full path to the file, without extension.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode write(const std::string &file) const = 0;
    
    /**
     * \brief read solution vectors from a given file.
     * 
     * Now only support HDF5 format.
     *
     * \param file [in] the full path to the file, without extension.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode read(const std::string &file) = 0;


protected:

    /**
     * \brief underlying initialization function.
     *
     * \param mesh [in] a type::Mesh.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode init(const type::Mesh &mesh) = 0;

    /** \brief MPI communicator. */
    MPI_Comm                                comm;

    /** \brief size of MPI communicator. */
    PetscMPIInt                             mpiSize;
    
    /** \brief the rank of this process. */
    PetscMPIInt                             mpiRank;

    /** \brief a shared_ptr to underlying mesh. */
    type::Mesh                              mesh;

};
} // end of namespace solution


namespace type
{
    /** \brief Solution type definition. */
    typedef std::shared_ptr<solution::SolutionBase> Solution;
}


namespace solution
{
    /**
     * \brief a factory function for creating Solution objects.
     *
     * \param mesh [in] a type::Mesh.
     * \param solution [out] resulting Solution object.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode createSolution(
            const type::Mesh &mesh, type::Solution &solution);
}

} // end of namespace petibm
