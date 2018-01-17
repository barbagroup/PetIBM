/**
 * \file boundary.h
 * \brief boundary::BoundaryBase, type::Boundary, and factory function.
 * \author Anush Krishnan (anus@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \copyright MIT.
 */


# pragma once

// STL
# include <memory>

// here goes headers from our PetIBM
# include <petibm/mesh.h>
# include <petibm/solution.h>
# include <petibm/singleboundary.h>


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
 * \f$a_0\f$ and \f$a_1\f$ are coefficients that depends on the types of boundary
 * conditions and also the location of the boundaries. For example, on a 
 * staggered grid, the relationship of a ghost point and its corresponding 
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
 * The data in a petibm::type::SingleBoundary are only distributed to the processes
 * owing ghost points on this boundary. And a member function of a 
 * petibm::type::SingleBoundary is called only on these processes. Except the
 * corners, different petibm::type::SingleBoundary are normally owned by different
 * process sets, so their member functions can be called concurrently.
 * 
 * While a petibm::type::SingleBoundary represents a single geometric boundary,
 * a petibm::type::Boundary represents a collection of all petibm::type::SingleBoundary
 * in a domain. For example, a petibm::type::Boundary of a 2D Cartesian mesh will
 * have 4 petibm::type::SingleBoundary in it. The design is, API users should use
 * the member functions in a petibm::type::Boundary to launch the functions of
 * all petibm::type::SingleBoundary in parallel. This can reduce the idle
 * time of CPU cores when doing something regarding boundaries.
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
 * \see boundaryModule, petibm::type::SingleBoundary, petibm::boundary::createBoundary
 * \ingroup boundaryModule
 */
class BoundaryBase
{
public:

    /** \brief Dimension. */
    PetscInt        dim;

    /** \brief A 2D vector holding all single boundaries. */
    std::vector<std::vector<type::SingleBoundary>>  bds;
    

    /** \brief Default constructor. */
    BoundaryBase() = default;

    /**
     * \brief Construct a boundary object based on a given mesh object.
     * \param mesh [in] a Mesh object.
     * \param node [in] a YAML::Node object.
     */
    BoundaryBase(const type::Mesh &mesh, const YAML::Node &node);

    /** \brief Default destructor. */
    virtual ~BoundaryBase();


    /**
     * \brief Manually destroy data.
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode destroy();


    /**
     * \brief Set the initial values of ghost points.
     * \param soln [in] a Solution object which will be used to calculate the values.
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode setGhostICs(const type::Solution &soln) = 0;

    /**
     * \brief Update the euqations between ghost and boundary points.
     * \param soln [in] a Solution object.
     * \param dt [in] time-step size.
     * \return PetscErrorCode.
     * 
     * Some kinds of boundary conditions will require changing equation coefficients.
     */
    virtual PetscErrorCode updateEqs(const type::Solution &soln, const PetscReal &dt) = 0;

    /**
     * \brief Update the values of ghost points.
     * \param soln [in] a Solution object.
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode updateGhostValues(const type::Solution &soln) = 0;

    /**
     * \brief Copy values of ghost points to a vector of local PETSc Vec objects.
     * \param lclVecs [in] a std::vector<Vec> object.
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode copyValues2LocalVecs(std::vector<Vec> &lclVecs) const = 0;

protected:

    /**
     * \brief Underlying initialization function.
     * \param mesh [in] a Mesh object.
     * \param node [in] a YAML::Node object.
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode init(const type::Mesh &mesh, const YAML::Node &node) = 0;

    /** \brief MPI communicator. */
    MPI_Comm                            comm;

    /** \brief Size of MPI communicator. */
    PetscMPIInt                         mpiSize;
    
    /** \brief The rank of this process. */
    PetscMPIInt                         mpiRank;

    /** \brief A shared_ptr to underlying mesh. */
    type::Mesh                          mesh;

};
} // end of namespace boundary


namespace type
{
    /**
     * \brief Type definition of petibm::type::Boundary.
     * \see boundaryModule, petibm::boundary::BoundaryBase, petibm::boundary::createBoundary
     * \ingroup boundaryModule
     * 
     * Please use petibm::boundary::createBoundary to create a Mesh object.
     */
    typedef std::shared_ptr<boundary::BoundaryBase> Boundary;
} // end of namespace type


namespace boundary
{
    /**
     * \brief Create a Boundary object.
     * \param mesh [in] underlying Mesh object.
     * \param node [in] a YAML::Node object.
     * \param boundary [out] resulting Boundary object.
     * \return PetscErrorCode.
     * \see boundaryModule, petibm::type::Boundary
     * \ingroup boundaryModule
     */
    PetscErrorCode createBoundary(const type::Mesh &mesh, 
            const YAML::Node &node, type::Boundary &boundary);
} // end of namespace boundary

} // end of namespace petibm
