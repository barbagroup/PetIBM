/**
 * \file boundary.h
 * \brief boundary::BoundaryBase, type::Boundary, and factory function.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 */

#pragma once

#include <memory>

#include <petibm/mesh.h>
#include <petibm/singleboundary.h>
#include <petibm/solution.h>

/**
 * \defgroup boundaryModule Boundaries and boundary conditions
 * \brief Objects handling boundary conditions with ghost point schemes.
 *
 * This module contains objects handling boundary conditions using ghost point
 * schemes. The concept is that each petibm::type::SingleBoundary holds the
 * information of the ghost points and corresponding boundary points on a
 * geometrical boundary. A petibm::type::SingleBoundary also holds the
 * information of the relationship between a ghost point and its corresponding
 * target point. The so called "relationship" means the coefficients of the
 * math equation connecting ghost points and the corresponding points. If we
 * denote the u velocity of a ghost point as \f$u_g\f$ and that of the
 * corresponding boundary point as \f$u_b\f$, then the relationship at time step
 * \f$t\f$ between these two can be written as
 * \f[u_g^t = a_0 \times u_b^t + a_1\f]
 * \f$a_0\f$ and \f$a_1\f$ are coefficients that depends on the types of
 * boundary conditions and also the location of the boundaries. For example, on
 * a staggered grid, the relationship of a ghost point and its corresponding
 * boundary point on \ref petibm::type::XPLUS "XPLUS" with Dirichlet
 * boundary condition at time step \f$t\f$ is
 * \f[u_g^t = 0 \times u_b^t + U_{Dirichlet}(t, y, z)\f]
 * Or for a convective boundary condition on \ref petibm::type::YPLUS
 * "YPLUS", the relationship will be
 * \f[
 * u_g^t = -1.0 \times u_b^t + \left(u_g^{t-1} + u_b^{t-1} -
 *      2.0 * \Delta t * U_c\frac{u_g^{t-1} - u_b^{t-1}}{\Delta L}\right)
 * \f]
 *
 * The data in a petibm::type::SingleBoundary are only distributed to the
 * processes owing ghost points on this boundary. And a member function of a
 * petibm::type::SingleBoundary is called only on these processes. Except the
 * corners, different petibm::type::SingleBoundary are normally owned by
 * different process sets, so their member functions can be called concurrently.
 *
 * While a petibm::type::SingleBoundary represents a single geometric boundary,
 * a petibm::type::Boundary represents a collection of all
 * petibm::type::SingleBoundary in a domain. For example, a
 * petibm::type::Boundary of a 2D Cartesian mesh will have 4
 * petibm::type::SingleBoundary in it. The design is, API users should use the
 * member functions in a petibm::type::Boundary to launch the functions of all
 * petibm::type::SingleBoundary in parallel. This can reduce the idle time of
 * CPU cores when doing something regarding boundaries.
 *
 *
 * \see petibm::type::Boundary, petibm::boundary::createBoundary
 * \ingroup petibm
 */

namespace petibm
{
/**
 * \brief A collection of all boundary-related objects and functions.
 * \see boundaryModule, petibm::type::Boundary, petibm::type::SingleBoundary
 * \ingroup boundaryModule
 */
namespace boundary
{
/**
 * \brief Base (abstract) class for the whole boundary.
 * \see boundaryModule, petibm::type::SingleBoundary,
 * petibm::boundary::createBoundary \ingroup boundaryModule
 */
class BoundaryBase
{
public:
    /** \brief Dimension. */
    PetscInt dim;

    /** \brief A 2D vector holding all single boundaries. */
    std::vector<std::vector<type::SingleBoundary>> bds;

    /** \brief Default constructor. */
    BoundaryBase() = default;

    /**
     * \brief Construct a boundary object based on a given mesh object.
     *
     * \param mesh [in] Structured Cartesian mesh object.
     * \param node [in] YAML configurations.
     */
    BoundaryBase(const type::Mesh &mesh, const YAML::Node &node);

    /** \brief Default destructor. */
    virtual ~BoundaryBase();

    /** \brief Manually destroy data.  */
    virtual PetscErrorCode destroy();

    /**
     * \brief Set the initial values of ghost points.
     *
     * \param soln [in] Solution object used to calculate the field solutions.
     */
    virtual PetscErrorCode setGhostICs(const type::Solution &soln) = 0;

    /**
     * \brief Update the equations between ghost and boundary points.
     *
     * \param soln [in] Data object for the Eulerian field solutions.
     * \param dt [in] Time-step size.
     *
     * Some kinds of boundary conditions will require changing equation
     * coefficients.
     */
    virtual PetscErrorCode updateEqs(const type::Solution &soln,
                                     const PetscReal &dt) = 0;

    /**
     * \brief Update the values of ghost points.
     *
     * \param soln [in] Data object for the Eulerian field solutions.
     */
    virtual PetscErrorCode updateGhostValues(const type::Solution &soln) = 0;

    /**
     * \brief Copy values of ghost points to a vector of local PETSc Vec
     * objects.
     *
     * \param lclVecs [in] a std::vector<Vec> object.
     */
    virtual PetscErrorCode copyValues2LocalVecs(
        std::vector<Vec> &lclVecs) const = 0;

protected:
    /**
     * \brief Underlying initialization function.
     *
     * \param mesh [in] Structured Cartesian mesh object.
     * \param node [in] YAML configurations.
     */
    virtual PetscErrorCode init(const type::Mesh &mesh,
                                const YAML::Node &node) = 0;

    /** \brief MPI communicator. */
    MPI_Comm comm;

    /** \brief Size of MPI communicator. */
    PetscMPIInt mpiSize;

    /** \brief The rank of this process. */
    PetscMPIInt mpiRank;

    /** \brief A shared_ptr to underlying mesh. */
    type::Mesh mesh;

};  // BoundaryBase

}  // end of namespace boundary

namespace type
{
/**
 * \brief Type definition of petibm::type::Boundary.
 *
 * \see boundaryModule, petibm::boundary::BoundaryBase,
 * petibm::boundary::createBoundary \ingroup boundaryModule
 *
 * Please use petibm::boundary::createBoundary to create a Mesh object.
 */
typedef std::shared_ptr<boundary::BoundaryBase> Boundary;

}  // end of namespace type

namespace boundary
{
/**
 * \brief Create a Boundary object.
 *
 * \param mesh [in] Structured Cartesian mesh object.
 * \param node [in] YAML configurations.
 * \param boundary [out] Data object with boundary conditions.
 *
 * \see boundaryModule, petibm::type::Boundary
 * \ingroup boundaryModule
 */
PetscErrorCode createBoundary(const type::Mesh &mesh, const YAML::Node &node,
                              type::Boundary &boundary);

}  // end of namespace boundary

}  // end of namespace petibm
