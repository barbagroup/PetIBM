/***************************************************************************//**
 * \file mesh.h
 * \author Anush Krishnan (anus@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Definition of the class `MeshBase` and factory function.
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


namespace petibm
{
namespace mesh
{
/** \brief base class of meshes. */
class MeshBase
{
public:

    /** \brief dimension. */
    PetscInt                dim = -1;

    /** \brief minimum coordinates of boundaries in all directions. */
    type::RealVec1D         min;

    /** \brief maximum coordinates of boundaries in all directions. */
    type::RealVec1D         max;

    /** \brief total number of points of all fields and in all directions. */
    type::IntVec2D          n;

    /** \brief coordinates of mesh points of all fields and in all directions. */
    type::GhostedVec3D      coord;

    /** \brief point spacing of mesh points of all fields and in all directions. */
    type::GhostedVec3D      dL;

    /** \brief total number of velocity points. */
    PetscInt                UN;

    /** \brief total number of pressure points. */
    PetscInt                pN;

    /** \brief a string for printing information. */
    std::string             info;



    // PETSc stuffs
    /** \brief a vector of DMs of all fields. */
    std::vector<DM>         da;

    /** \brief number of processes in all directions. */
    type::IntVec1D          nProc;

    /** \brief the beginging index of all fields in all directions of this process. */
    type::IntVec2D          bg;

    /** \brief the ending index of all fields in all directions of this process. */
    type::IntVec2D          ed;

    /** \brief the number of points of all fields in all directions of this process. */
    type::IntVec2D          m;

    /** \brief total number of velocity points local to this process. */
    PetscInt                UNLocal;

    /** \brief total number of pressure points local to this process. */
    PetscInt                pNLocal;

    /** \brief DMComposte of velocity DMs. */
    DM                      UPack;


    // MPI stuffs
    /** \brief communicator. */
    MPI_Comm                                comm;

    /** \brief total number of processes. */
    PetscMPIInt                             mpiSize;

    /** \brief rank of this process. */
    PetscMPIInt                             mpiRank;



    /** \brief default constructor. */
    MeshBase() = default;


    /**
     * \brief constructor.
     *
     * \param world MPI communicator.
     * \param node a YAML node containing setting of Cartesian mesh.
     */
    MeshBase(const MPI_Comm &world, const YAML::Node &node) {};


    /** \brief default destructor. */
    virtual ~MeshBase() = default;
    

    /**
     * \brief print iformation to standard output.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode printInfo() const;


    /**
     * \brief get the local index of a point by providing MatStencil.
     *
     * Note that the stencil must belong to this process.
     *
     * \param f target field (u=0, v=1, w=2, p=3).
     * \param s MatStencil of target velocity point.
     * \param idx returned local index.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode getLocalIndex(
            const PetscInt &f, const MatStencil &s, PetscInt &idx) const = 0;


    /**
     * \brief get the local index of a point by i, j, k.
     *
     * Note that the stencil must belong to this process. For 2D mesh, the
     * value of k-index will be ignored.
     *
     * \param f target field (u=0, v=1, w=2, p=3).
     * \param i i-index.
     * \param j j-index.
     * \param k k-index.
     * \param idx returned local index.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode getLocalIndex(const PetscInt &f, 
            const PetscInt &i, const PetscInt &j, const PetscInt &k, 
            PetscInt &idx) const = 0;


    /**
     * \brief get the natural index of a point by providing MatStencil.
     *
     * \param f target field (u=0, v=1, w=2, p=3).
     * \param s MatStencil of target velocity point.
     * \param idx returned natural index.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode getNaturalIndex(
            const PetscInt &f, const MatStencil &s, PetscInt &idx) const = 0;


    /**
     * \brief get the natural index of a point by i, j, k.
     *
     * For 2D mesh, the value of k-index will be ignored.
     *
     * \param f target field (u=0, v=1, w=2, p=3).
     * \param i i-index.
     * \param j j-index.
     * \param k k-index.
     * \param idx returned natural index.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode getNaturalIndex(const PetscInt &f, 
            const PetscInt &i, const PetscInt &j, const PetscInt &k, 
            PetscInt &idx) const = 0;


    /**
     * \brief get the global index of a point in un-packed DM by 
     *        providing MatStencil.
     *
     * \param f target field (u=0, v=1, w=2, p=3).
     * \param s MatStencil of target velocity point.
     * \param idx returned global index in un-packed DM.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode getGlobalIndex(
            const PetscInt &f, const MatStencil &s, PetscInt &idx) const = 0;


    /**
     * \brief get the global index of a point in un-packed DM by 
     *        i, j, k.
     *
     * For 2D mesh, the value of k-index will be ignored.
     *
     * \param f target field (u=0, v=1, w=2, p=3).
     * \param i i-index.
     * \param j j-index.
     * \param k k-index.
     * \param idx returned global index in un-packed DM.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode getGlobalIndex(const PetscInt &f, 
            const PetscInt &i, const PetscInt &j, const PetscInt &k, 
            PetscInt &idx) const = 0;


    /**
     * \brief get the global index of a point in packed DM by providing MatStencil.
     * 
     * For pressure field (f=3), packed ID will be equal to unpacked ID.
     *
     * \param f target field (u=0, v=1, w=2, p=3).
     * \param s MatStencil of target velocity point.
     * \param idx returned global index in packed DM.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode getPackedGlobalIndex(
            const PetscInt &f, const MatStencil &s, PetscInt &idx) const = 0;


    /**
     * \brief get the global index of a point in packed DM by 
     *        i, j, k.
     *
     * For 2D mesh, the value of k-index will be ignored. For pressure field 
     * (f=3), packed ID will be equal to unpacked ID.
     *
     * \param f target field (u=0, v=1, w=2, p=3).
     * \param i i-index.
     * \param j j-index.
     * \param k k-index.
     * \param idx returned global index in packed DM.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode getPackedGlobalIndex(const PetscInt &f, 
            const PetscInt &i, const PetscInt &j, const PetscInt &k, 
            PetscInt &idx) const = 0;

protected:


    /**
     * \brief initialization.
     *
     * \param world MPI communicator.
     * \param node a YAML node containing setting of Cartesian mesh.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode init(
            const MPI_Comm &world, const YAML::Node &node) = 0;
    
}; // MeshBase
} // end of namespace mesh


namespace type
{
    /** \brief Mesh type definition. */
    typedef std::shared_ptr<mesh::MeshBase> Mesh;
} // end of type


namespace mesh
{
    /**
     * \brief factory function for creating a Mesh object.
     *
     * \param node [in] YAML node.
     *
     * \return Mesh object (a shared_ptr to MeshBase).
     */
    PetscErrorCode createMesh(
            const MPI_Comm &comm, const YAML::Node &node, type::Mesh &mesh);
} // end of namespace mesh

} // end of namespace petibm
