/***************************************************************************//**
 * \file boundary.h
 * \author Anush Krishnan (anus@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Definition of the abstract class `BoundaryBase` and related utilities.
 */


# pragma once

// STL
# include <memory>

// here goes headers from our PetIBM
# include <petibm/mesh.h>
# include <petibm/solution.h>


namespace petibm
{
namespace boundary
{
/** \brief base abstract class for the whole boundary. */
class BoundaryBase
{
public:

    /** \brief dimension. */
    PetscInt        dim;
    

    /** \brief default constructor. */
    BoundaryBase() = default;

    /**
     * \brief construct a boundary object based on a given mesh object.
     *
     * \param mesh [in] a Mesh object.
     * \param node [in] a YAML::Node object.
     */
    BoundaryBase(const type::Mesh &mesh, const YAML::Node &node);

    /** \brief default destructor. */
    virtual ~BoundaryBase() = default;


    /**
     * \brief set the initial values of ghost points.
     *
     * \param soln [in] a Solution object which will be used to calculate the values.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode setGhostICs(const type::Solution &soln) = 0;

    /**
     * \brief update the euqations between ghost and boundary points.
     * 
     * Some kinds of boundary conditions will require changing equation coefficients.
     *
     * \param soln [in] a Solution object.
     * \param dt [in] time-step size.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode updateEqs(const type::Solution &soln, const PetscReal &dt) = 0;

    /**
     * \brief update the values of ghost points.
     *
     * \param soln [in] a Solution object.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode updateGhostValues(const type::Solution &soln) = 0;

    /**
     * \brief copy values of ghost points to a vector of local PETSc Vec objects.
     *
     * \param lclVecs [in] a std::vector<Vec> object.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode copyValues2LocalVecs(std::vector<Vec> &lclVecs) const = 0;

protected:

    /**
     * \brief underlying initialization function.
     *
     * \param mesh [in] a Mesh object.
     * \param node [in] a YAML::Node object.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode init(const type::Mesh &mesh, const YAML::Node &node) = 0;

    /** \brief MPI communicator. */
    MPI_Comm                            comm;

    /** \brief size of MPI communicator. */
    PetscMPIInt                         mpiSize;
    
    /** \brief the rank of this process. */
    PetscMPIInt                         mpiRank;

    /** \brief a shared_ptr to underlying mesh. */
    type::Mesh                          mesh;

};
} // end of namespace boundary


namespace type
{
    /** \brief Boundary type definition. */
    typedef std::shared_ptr<boundary::BoundaryBase> Boundary;
} // end of namespace type


namespace boundary
{
    /**
     * \brief create a Boundary object.
     *
     * \param mesh [in] underlying Mesh object.
     * \param node [in] a YAML::Node object.
     * \param boundary [out] resulting Boundary object.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode createBoundary(const type::Mesh &mesh, 
            const YAML::Node &node, type::Boundary &boundary);
} // end of namespace boundary

} // end of namespace petibm
