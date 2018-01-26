/**
 * \file mesh.h
 * \brief Prototype of mesh::MeshBase, type::Mesh, and factory function.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 */


# pragma once

// here goes C++ STL
# include <string>
# include <vector>

// here goes PETSc headers
# include <petscsys.h>
# include <petscao.h>
# include <petscdmda.h>
# include <petscdmcomposite.h>

// here goes YAML header
# include <yaml-cpp/yaml.h>

// here goes headers from our PetIBM
# include <petibm/type.h>


/**
 * \defgroup meshModule Mesh objects
 * \brief Objects for meshes and domain decompositions
 * 
 * So far, we only have stretched Cartesian mesh in our code. However, in order
 * to consider future extension of different kinds of mesh, we still design a
 * base (abstract) class, \ref petibm::mesh::MeshBase "MeshBase", to unify the 
 * interfaces of potential mesh classes. Currently, there is only a child 
 * derived from the base class: \ref petibm::mesh::CartesianMesh "CartesianMesh".
 * 
 * The design of PetIBM is to use `std::shared_ptr` to hold the base class for
 * instances of different kinds of derived classes, so that we can call the public 
 * APIs regardless the types of the derived classes. \ref petibm::type::Mesh "Mesh" 
 * is defined as a shared pointer to \ref petibm::mesh::MeshBase "MeshBase". Users
 * should use \ref petibm::type::Mesh "Mesh" instead of using classes directly. 
 * And users should also initialize mesh instances with the factory function, 
 * \ref petibm::mesh::createMesh "createMesh", instead of initializing the 
 * instances directly.
 * 
 * The domain decomposition is done with PETSc DMDA objects. And the information
 * of sub-domains are retrieved into mesh classes.
 * 
 * \see petibm::type::Mesh, petibm::mesh::createMesh
 * \ingroup petibm
 */


namespace petibm
{
/** 
 * \brief Collection of classes and utilities regarding mesh.
 * \see meshModule, petibm::type::Mesh, petibm::mesh::createMesh
 * \ingroup meshModule
 */
namespace mesh
{
/**
 * \brief Base (abstract) class of mesh.
 * \see meshModule, petibm::type::Mesh, petibm::mesh::createMesh
 * \ingroup meshModule
 */
class MeshBase
{
public:

    /** \brief Dimension. */
    PetscInt                dim = -1;

    /** \brief Minimum coordinates of boundaries in all directions. */
    type::RealVec1D         min;

    /** \brief Maximum coordinates of boundaries in all directions. */
    type::RealVec1D         max;

    /** \brief Total number of points of all fields and in all directions. */
    type::IntVec2D          n;
    
    /** \brief Bools indicating if any direction is periodic. */
    type::BoolVec2D         periodic;

    /** \brief Coordinates of mesh points of all fields and in all directions. */
    type::GhostedVec3D      coord;

    /** \brief Spacings of mesh points of all fields and in all directions. */
    type::GhostedVec3D      dL;

    /** \brief Total number of velocity points. */
    PetscInt                UN;

    /** \brief Total number of pressure points. */
    PetscInt                pN;

    /** \brief A string for printing information. */
    std::string             info;



    // PETSc stuffs
    /** \brief A vector of DMs of all fields. */
    std::vector<DM>         da;

    /** \brief Number of processes in all directions. */
    type::IntVec1D          nProc;

    /** \brief The beginning index of all fields in all directions of this process. */
    type::IntVec2D          bg;

    /** \brief The ending index of all fields in all directions of this process. */
    type::IntVec2D          ed;

    /** \brief The number of points of all fields in all directions of this process. */
    type::IntVec2D          m;

    /** \brief Total number of velocity points local to this process. */
    PetscInt                UNLocal;

    /** \brief Total number of pressure points local to this process. */
    PetscInt                pNLocal;

    /** \brief DMComposte of velocity DMs. */
    DM                      UPack;


    // MPI stuffs
    /** \brief Communicator. */
    MPI_Comm                                comm;

    /** \brief Total number of processes. */
    PetscMPIInt                             mpiSize;

    /** \brief Rank of this process. */
    PetscMPIInt                             mpiRank;



    /** \brief Default constructor. */
    MeshBase() = default;


    /**
     * \brief Constructor.
     *
     * \param world [in] MPI communicator.
     * \param node [in] a YAML::node containing settings of mesh.
     * 
     * Users are not encouraged to use constructor to initialize instances
     * directly. Please use the factory function.
     * 
     * \see petibm::mesh::createMesh
     */
    MeshBase(const MPI_Comm &world, const YAML::Node &node) {};


    /** \brief Default destructor. */
    virtual ~MeshBase();


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
     * \brief Write grid/mesh data to a HDF5.
     *
     * \param filePath [in] path to the file (without extension).
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode write(const std::string &filePath) const = 0;


    /**
     * \brief Get the local index of a point by providing MatStencil.
     *
     * \param f [in] target field (u=0, v=1, w=2, p=3).
     * \param s [in] MatStencil of target point.
     * \param idx [out] local index.
     *
     * Note that the stencil must belong to this process.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode getLocalIndex(
            const PetscInt &f, const MatStencil &s, PetscInt &idx) const = 0;


    /**
     * \brief Get the local index of a point by providing i, j, and k.
     *
     * \param f [in] target field (u=0, v=1, w=2, p=3).
     * \param i [in] i-index.
     * \param j [in] j-index.
     * \param k [in] k-index.
     * \param idx [out] local index.
     *
     * Note that the stencil must belong to this process. For 2D mesh, the
     * value of k-index will be ignored.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode getLocalIndex(const PetscInt &f, 
            const PetscInt &i, const PetscInt &j, const PetscInt &k, 
            PetscInt &idx) const = 0;


    /**
     * \brief Get the natural index of a point by providing MatStencil.
     *
     * \param f [in] target field (u=0, v=1, w=2, p=3).
     * \param s [in] MatStencil of target point.
     * \param idx [out] returned natural index.
     * 
     * If the provided MatStencil is not valid or is a ghost point, the `idx`
     * will be `-1`.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode getNaturalIndex(
            const PetscInt &f, const MatStencil &s, PetscInt &idx) const = 0;


    /**
     * \brief Get the natural index of a point by providing i, j, and k.
     *
     * \param f [in] target field (u=0, v=1, w=2, p=3).
     * \param i [in] i-index.
     * \param j [in] j-index.
     * \param k [in] k-index.
     * \param idx [out] natural index.
     *
     * For 2D mesh, the value of k-index will be ignored.
     * 
     * If the provided `(i, j, k)` is not valid or is a ghost point, the `idx`
     * will be `-1`.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode getNaturalIndex(const PetscInt &f, 
            const PetscInt &i, const PetscInt &j, const PetscInt &k, 
            PetscInt &idx) const = 0;


    /**
     * \brief Get the global index of a point in unpacked DM by 
     *        providing MatStencil.
     *
     * \param f [in] target field (u=0, v=1, w=2, p=3).
     * \param s [in] MatStencil of target point.
     * \param idx [out] global index in unpacked DM.
     * 
     * If the provided MatStencil is not valid or is a ghost point, the `idx`
     * will be `-1`.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode getGlobalIndex(
            const PetscInt &f, const MatStencil &s, PetscInt &idx) const = 0;


    /**
     * \brief Get the global index of a point in unpacked DM by providing
     *        i, j, and k.
     *
     * \param f [in] target field (u=0, v=1, w=2, p=3).
     * \param i [in] i-index.
     * \param j [in] j-index.
     * \param k [in] k-index.
     * \param idx [out] global index in unpacked DM.
     *
     * For 2D mesh, the value of k-index will be ignored.
     * 
     * If the provided `(i, j, k)` is not valid or is a ghost point, the `idx`
     * will be `-1`.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode getGlobalIndex(const PetscInt &f, 
            const PetscInt &i, const PetscInt &j, const PetscInt &k, 
            PetscInt &idx) const = 0;


    /**
     * \brief Get the global index of a point in packed DM by providing MatStencil.
     *
     * \param f [in] target field (u=0, v=1, w=2, p=3).
     * \param s [in] MatStencil of target point.
     * \param idx [out] global index in packed DM.
     * 
     * For pressure field (f=3), packed ID will be equal to unpacked ID.
     * 
     * If the provided MatStencil is not valid or is a ghost point, the `idx`
     * will be `-1`.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode getPackedGlobalIndex(
            const PetscInt &f, const MatStencil &s, PetscInt &idx) const = 0;


    /**
     * \brief Get the global index of a point in packed DM by providing
     *        i, j, and k.
     *
     * \param f [in] target field (u=0, v=1, w=2, p=3).
     * \param i [in] i-index.
     * \param j [in] j-index.
     * \param k [in] k-index.
     * \param idx [out] global index in packed DM.
     *
     * For 2D mesh, the value of k-index will be ignored. For pressure field 
     * (f=3), packed ID will be equal to unpacked ID.
     * 
     * If the provided `(i, j, k)` is not valid or is a ghost point, the `idx`
     * will be `-1`.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode getPackedGlobalIndex(const PetscInt &f, 
            const PetscInt &i, const PetscInt &j, const PetscInt &k, 
            PetscInt &idx) const = 0;

protected:


    /**
     * \brief Initialization.
     *
     * \param world [in] MPI communicator.
     * \param node [in] a YAML node containing setting of Cartesian mesh.
     * 
     * Given the design of PetIBM is to use factory function to create shared
     * pointers of instances, it seems not necessary to make this function
     * public.
     *
     * \return PetscErrorCode.
     * 
     * \see petibm::mesh::createMesh
     */
    virtual PetscErrorCode init(
            const MPI_Comm &world, const YAML::Node &node) = 0;
    
}; // MeshBase
} // end of namespace mesh


namespace type
{
    /**
     * \brief Type definition of Mesh.
     * 
     * Please use petibm::mesh::createMesh to create a Mesh object.
     * 
     * Example usage:
     * \code
     * PetscErrorCode ierr;
     * petibm::type::Mesh mesh;
     * 
     * // create Mesh with petibm::mesh::createMesh
     * 
     * ierr = mesh->printInfo(); CHKERRQ(ierr);
     * ierr = mesh->write("./grid"); CHKERRQ(ierr);
     * 
     * PetscInt idx;
     * ierr = mesh->getNaturalIndex(0, 1, 1, 1, idx); CHKERRQ(ierr);
     * \endcode
     * 
     * \see meshModule, petibm::mesh::createMesh, petibm::mesh::MeshBase
     * \ingroup meshModule
     */
    typedef std::shared_ptr<mesh::MeshBase> Mesh;
} // end of namespace type


namespace mesh
{
    /**
     * \brief Factory function for creating a Mesh object.
     * \param comm [in] MPI Communicator.
     * \param node [in] a YAML node containing all settings.
     * \param mesh [out] a type::Mesh instance.
     * \return PetscErrorCode.
     * 
     * This function will look into the key `mesh` in `node` for mesh settings. 
     * An example of 2D stretched Cartesian mesh settings in the YAML node is 
     * shown as the following:
     * 
     * \code
     * YAML::Node node;
     * node["mesh"][0]["direction"] = "x";
     * node["mesh"][0]["start"] = 0.0;
     * node["mesh"][0]["subDomains"][0]["end"] = 0.5;
     * node["mesh"][0]["subDomains"][0]["cells"] = 64;
     * node["mesh"][0]["subDomains"][0]["stretchRatio"] = 1.01;
     * node["mesh"][0]["subDomains"][1]["end"] = 1.0;
     * node["mesh"][0]["subDomains"][1]["cells"] = 64;
     * node["mesh"][0]["subDomains"][1]["stretchRatio"] = 0.9900990099001;
     * node["mesh"][1]["direction"] = "y";
     * node["mesh"][1]["start"] = 0.0;
     * node["mesh"][1]["subDomains"][0]["end"] = 0.5;
     * node["mesh"][1]["subDomains"][0]["cells"] = 64;
     * node["mesh"][1]["subDomains"][0]["stretchRatio"] = 1.01;
     * node["mesh"][1]["subDomains"][1]["end"] = 1.0;
     * node["mesh"][1]["subDomains"][1]["cells"] = 64;
     * node["mesh"][1]["subDomains"][1]["stretchRatio"] = 0.9900990099001;
     * \endcode
     * 
     * Or in an ASCII YAML file:
     * 
     * ```
     * mesh:
     *   - direction: x
     *     start: 0.0
     *     subDomains:
     *       - end: 0.5
     *         cells: 64
     *         stretchRatio: 1.01
     *       - end: 1.0
     *         cells: 64
     *         stretchRatio: 0.9900990099001
     *   
     *   - direction: y
     *     start: 0.0
     *     subDomains:
     *       - end: 0.5
     *         cells: 64
     *         stretchRatio: 1.01
     *       - end: 1.0
     *         cells: 64
     *         stretchRatio: 0.9900990099001
     * ```
     * 
     * We then pass this YAML node into the factory function:
     * \code
     * PetscErrorCode ierr;
     * petibm::type::mesh mesh;
     * ierr = petibm::mesh::createMesh(PETSC_COMM_WORLD, node, mesh); CHKERRQ(ierr);
     * \endcode
     * 
     * \see meshModule, petibm::type::Mesh
     * \ingroup meshModule
     */
    PetscErrorCode createMesh(
            const MPI_Comm &comm, const YAML::Node &node, type::Mesh &mesh);
} // end of namespace mesh

} // end of namespace petibm
