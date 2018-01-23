/**
 * \file cartesianmesh.h
 * \brief Definition of class mesh::CartesianMesh.
 * \author Anush Krishnan (anush@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \copyright MIT.
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
# include <petibm/mesh.h>


namespace petibm
{
namespace mesh
{

/**
 * \brief Class of composite staggered Cartesian mesh.
 * \see meshModule, petibm::type::Mesh, petibm::mesh::createMesh
 * \ingroup meshModule
 */
class CartesianMesh : public MeshBase
{
public:

    /** \copydoc MeshBase(const MPI_Comm &, const YAML::Node &) */
    CartesianMesh(const MPI_Comm &world, const YAML::Node &node);

    /** \copydoc ~MeshBase */
    virtual ~CartesianMesh();


    // doc is the same as MeshBase::destroy
    virtual PetscErrorCode destroy();
    
    // doc is the same as MeshBase::write
    virtual PetscErrorCode write(const std::string &filePath) const;

    // doc is the same as MeshBase::getLocalIndex
    virtual PetscErrorCode getLocalIndex(
            const PetscInt &f, const MatStencil &s, PetscInt &idx) const;

    // doc is the same as MeshBase::getLocalIndex
    virtual PetscErrorCode getLocalIndex(const PetscInt &f, 
            const PetscInt &i, const PetscInt &j, const PetscInt &k, 
            PetscInt &idx) const;

    // doc is the same as MeshBase::getNaturalIndex
    virtual PetscErrorCode getNaturalIndex(
            const PetscInt &f, const MatStencil &s, PetscInt &idx) const;

    // doc is the same as MeshBase::getNaturalIndex
    virtual PetscErrorCode getNaturalIndex(const PetscInt &f, 
            const PetscInt &i, const PetscInt &j, const PetscInt &k, 
            PetscInt &idx) const;

    // doc is the same as MeshBase::getGlobalIndex
    virtual PetscErrorCode getGlobalIndex(
            const PetscInt &f, const MatStencil &s, PetscInt &idx) const;

    // doc is the same as MeshBase::getGlobalIndex
    virtual PetscErrorCode getGlobalIndex(const PetscInt &f, 
            const PetscInt &i, const PetscInt &j, const PetscInt &k, 
            PetscInt &idx) const;

    // doc is the same as MeshBase::getPackedGlobalIndex
    virtual PetscErrorCode getPackedGlobalIndex(
            const PetscInt &f, const MatStencil &s, PetscInt &idx) const;

    // doc is the same as MeshBase::getPackedGlobalIndex
    virtual PetscErrorCode getPackedGlobalIndex(const PetscInt &f, 
            const PetscInt &i, const PetscInt &j, const PetscInt &k, 
            PetscInt &idx) const;

protected:
    
    // doc is the same as MeshBase::init
    PetscErrorCode init(const MPI_Comm &world, const YAML::Node &node);


    /** \brief Underlying data for mesh point spacing. */
    type::RealVec3D         dLTrue;

    /** \brief Underlying data for mesh point coordinates. */
    type::RealVec3D         coordTrue;

    /** \brief References to underlying AO objects of velocity DMs. */
    std::vector<AO>         ao;

    /** \brief Number of local velocity points (without ghost points) for all 
     *         processes and all velocity fields. */
    type::IntVec2D          UNLocalAllProcs;

    /** \brief Number of local packed velocity points (without ghost points) for 
     *         all processes. */
    type::IntVec1D          UPackNLocalAllProcs;

    /** \brief Offsets of velocity points in unpacked DMs for all processes and
     *         all velocity field. */
    type::IntVec2D          offsetsAllProcs;

    /** \brief Offsets of packed velocity points in packed DM. */
    type::IntVec1D          offsetsPackAllProcs;

    
    /**
     * \brief Create vertex information.
     * \return PetscErrorCode.
     */
    PetscErrorCode createVertexMesh();

    /**
     * \brief Create pressure mesh information.
     * \return PetscErrorCode.
     */
    PetscErrorCode createPressureMesh();

    /**
     * \brief Create velocity mesh information.
     * \return PetscErrorCode.
     */
    PetscErrorCode createVelocityMesh();

    /**
     * \brief Create a string of information.
     * \return PetscErrorCode.
     */
    PetscErrorCode createInfoString();

    /**
     * \brief Gather information of parallel distribution and add to info string.
     * \return PetscErrorCode.
     */
    PetscErrorCode addLocalInfoString(std::stringstream &ss);

    /**
     * \brief Initialize DMDAs.
     * \return PetscErrorCode.
     */
    PetscErrorCode initDMDA();

    /**
     * \brief Function for creating a single DMDA.
     * \param i [in] The index of the targeting field (0 ~ 4 represents u, v, w, 
     *          pressure, and vertex respectively).
     * \return PetscErrorCode.
     */
    PetscErrorCode createSingleDMDA(const PetscInt &i);

    /**
     * \brief Create DMDA for pressure.
     * \return PetscErrorCode.
     */
    PetscErrorCode createPressureDMDA();

    /**
     * \brief Create DMDAs for velocity fields and make a DMComposite.
     * \return PetscErrorCode.
     */
    PetscErrorCode createVelocityPack();

}; // CartesianMesh

} // end of namespace mesh
} // end of namespace petibm
