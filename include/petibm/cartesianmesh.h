/***************************************************************************//**
 * \file CartesianMesh.h
 * \author Anush Krishnan (anus@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Definition of the class `CartesianMesh`.
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

/** \brief class of composite Cartesian meshes of different field. */
class CartesianMesh : public MeshBase
{
public:

    /** \copydoc petibm::mesh::MeshBase::MeshBase */
    CartesianMesh(const MPI_Comm &world, const YAML::Node &node);

    /** \brief default destructor. */
    virtual ~CartesianMesh();
    
    /** \copydoc petibm::mesh::MeshBase::write */
    virtual PetscErrorCode write(const std::string &filePath) const;

    /** \copydoc petibm::mesh::MeshBase::getLocalIndex */
    virtual PetscErrorCode getLocalIndex(
            const PetscInt &f, const MatStencil &s, PetscInt &idx) const;

    /** \copydoc petibm::mesh::MeshBase::getLocalIndex */
    virtual PetscErrorCode getLocalIndex(const PetscInt &f, 
            const PetscInt &i, const PetscInt &j, const PetscInt &k, 
            PetscInt &idx) const;

    /** \copydoc petibm::mesh::MeshBase::getNaturalIndex */
    virtual PetscErrorCode getNaturalIndex(
            const PetscInt &f, const MatStencil &s, PetscInt &idx) const;

    /** \copydoc petibm::mesh::MeshBase::getNaturalIndex */
    virtual PetscErrorCode getNaturalIndex(const PetscInt &f, 
            const PetscInt &i, const PetscInt &j, const PetscInt &k, 
            PetscInt &idx) const;

    /** \copydoc petibm::mesh::MeshBase::getGlobalIndex */
    virtual PetscErrorCode getGlobalIndex(
            const PetscInt &f, const MatStencil &s, PetscInt &idx) const;

    /** \copydoc petibm::mesh::MeshBase::getGlobalIndex */
    virtual PetscErrorCode getGlobalIndex(const PetscInt &f, 
            const PetscInt &i, const PetscInt &j, const PetscInt &k, 
            PetscInt &idx) const;

    /** \copydoc petibm::mesh::MeshBase::getPackedGlobalIndex */
    virtual PetscErrorCode getPackedGlobalIndex(
            const PetscInt &f, const MatStencil &s, PetscInt &idx) const;

    /** \copydoc petibm::mesh::MeshBase::getPackedGlobalIndex */
    virtual PetscErrorCode getPackedGlobalIndex(const PetscInt &f, 
            const PetscInt &i, const PetscInt &j, const PetscInt &k, 
            PetscInt &idx) const;

protected:
    
    /** \copydoc petibm::mesh::MeshBase::init */
    PetscErrorCode init(const MPI_Comm &world, const YAML::Node &node);


    /** \brief the underlying data for mesh point spacing. */
    type::RealVec3D         dLTrue;

    /** \brief the underlying data for mesh point coordinates. */
    type::RealVec3D         coordTrue;

    /** \brief references to underlying AO objects of velocity DMs. */
    std::vector<AO>         ao;

    /** \brief number of local velocity points (without ghosts) for all 
     *         processes and all velocity fields. */
    type::IntVec2D          UNLocalAllProcs;

    /** \brief number of local packed velocity points (without ghost) for 
     *         all processes. */
    type::IntVec1D          UPackNLocalAllProcs;

    /** \brief offsets of velocity points in un-packed DMs for each field. */
    type::IntVec2D          offsetsAllProcs;

    /** \brief offsets of packed velocity points in packed DM. */
    type::IntVec1D          offsetsPackAllProcs;

    
    /**
     * \brief create vertex information.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode createVertexMesh();


    /**
     * \brief create pressure mesh information.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode createPressureMesh();


    /**
     * \brief create velocity mesh information.
     * 
     * \param periodic [in] bools indicating periodicity of boundaries.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode createVelocityMesh();


    /**
     * \brief create a string of information.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode createInfoString();


    /**
     * \brief gather information of parallel distribution and add to info string.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode addLocalInfoString(std::stringstream &ss);


    /**
     * \brief initialize DMDAs.
     *
     * \param periodic [in] bools indicating periodicity of boundaries.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode initDMDA();


    /**
     * \brief function for creating a single DMDA.
     *
     * \param i the index of the targeting field (0 ~ 4 represents u, v, w, 
     *          pressure, and vertex respectively.
     * \param periodic [in] bools indicating periodicity of boundaries.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode createSingleDMDA(const PetscInt &i);


    /**
     * \brief create DMDA for pressure.
     *
     * \param periodic [in] bools indicating periodicity of boundaries.
     * 
     * \return PetscErrorCode.
     */
    PetscErrorCode createPressureDMDA();


    /**
     * \brief create DMDAs for velovity fields and make a DMComposite.
     *
     * \param periodic [in] bools indicating periodicity of boundaries.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode createVelocityPack();

}; // CartesianMesh

} // end of namespace mesh
} // end of namespace petibm
