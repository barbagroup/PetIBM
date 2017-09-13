/***************************************************************************//**
 * \file bodypack.h
 * \author Anush Krishnan (anus@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Definition of the class `BodyPack`.
 */


# pragma once


// STL
# include <string>
# include <vector>
# include <memory>

// PETSc
# include <petscsys.h>
# include <petscdm.h>
# include <petscis.h>

// YAML
# include <yaml-cpp/yaml.h>

// PetIBM
# include "singlebody.h"


namespace petibm
{
namespace body
{

// TODO: should we really need ISLocalToGlobalMapping?
/** \brief class for a pack of multiple bodies. */
class BodyPack
{
public:

    /** \brief dimension. */
    PetscInt                    dim;

    /** \brief number of bodies in this pack. */
    PetscInt                    nBodies;

    /** \brief total number of Lagrangian points. */
    PetscInt                    nPts;

    /** \brief total number of local Lagrangian points. */
    PetscInt                    nLclPts;

    /** \brief a vector of SingleBody instances. */
    std::vector<SingleBody>     bodies;

    /** \brief a DMComposite of DMs of all `SingleBody`s. */
    DM                          dmPack;

    /** \brief local to global (in DMComposite) mappings for bodies. */
    std::vector<ISLocalToGlobalMapping>     mapping;

    /** \brief a string for printing information. */
    std::string                 info;


    /** \brief default constructor. */
    BodyPack();


    /**
     * \brief constructor of using CartesianMesh and a YAML node.
     *
     * \param mesh an instance of CartesianMesh for background mesh.
     * \param node a YAML node specifying information of all bodies.
     */
    BodyPack(const CartesianMesh &mesh, const YAML::Node &node);


    /** \brief default destructor. */
    ~BodyPack();


    /**
     * \brief initialization using CartesianMesh and a YAML node.
     *
     * \param mesh an instance of CartesianMesh for background mesh.
     * \param node a YAML node specifying information of all bodies.
     */
    PetscErrorCode init(const CartesianMesh &mesh, const YAML::Node &node);


    /**
     * \brief print information.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode printInfo() const;


    /**
     * \brief find which process owns the target Lagrangian point of target body.
     *
     * \param bIdx index of target body.
     * \param ptIdx index of target point.
     * \param proc returned process index.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode findProc(const PetscInt &bIdx, 
            const PetscInt &ptIdx, PetscMPIInt &proc) const;


    /**
     * \brief find un-packed global index of a DoF of Lagrangian point of a body.
     *
     * \param bIdx index of target body.
     * \param ptIdx index of target point.
     * \param dof index of target DoF.
     * \param idx returned un-packed global index.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode getGlobalIndex(const PetscInt &bIdx, 
            const PetscInt &ptIdx, const PetscInt &dof, PetscInt &idx) const;


    /**
     * \brief find un-packed global index of a DoF of Lagrangian point of a body.
     *
     * \param bIdx index of target body.
     * \param s MatStencil of target point.
     * \param idx returned un-packed global index.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode getGlobalIndex(const PetscInt &bIdx, 
            const MatStencil &s, PetscInt &idx) const;



    /**
     * \brief find packed global index of a DoF of Lagrangian point of a body.
     *
     * \param bIdx index of target body.
     * \param ptIdx index of target point.
     * \param dof index of target DoF.
     * \param idx returned packed global index.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode getPackedGlobalIndex(const PetscInt &bIdx, 
            const PetscInt &ptIdx, const PetscInt &dof, PetscInt &idx) const;



    /**
     * \brief find packed global index of a DoF of Lagrangian point of a body.
     *
     * \param bIdx index of target body.
     * \param s MatStencil of target point.
     * \param idx returned packed global index.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode getPackedGlobalIndex(const PetscInt &bIdx, 
            const MatStencil &s, PetscInt &idx) const;


    /**
     * \brief calculate the averaged force of each body.
     *
     * \param f packed force Vec of Lagrangian points.
     * \param fAvg return averaged force for each body.
     *
     * Note: fAvg doesn't have to have correct size. This function will resizing
     * the fAvg.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode calculateAvgForces(const Vec &f, types::RealVec2D &fAvg);


    /**
     * \brief calculate averaged forces for each body and directly write to a file.
     *
     * \param t current simulation time.
     * \param f packed Lagrangian force Vec.
     * \param dir directory in which the file will be.
     * \param file the file name without extension.
     *
     * Note: the forces will be appended to the file if the file already exists.
     *
     * \return PetscErrorCode. 
     */
    PetscErrorCode writeAvgForce(const PetscReal &dt, const Vec &f, 
            const std::string &dir, const std::string &file);

protected:

    /** \brief reference to backgrounf CartesianMesh. */
    std::shared_ptr<const CartesianMesh>    mesh;


    /** \brief reference to the MPI communicator. */
    std::shared_ptr<const MPI_Comm>         comm;


    /** \brief the total number of processes. */
    PetscMPIInt                             mpiSize;


    /** \brief the rank of this process. */
    PetscMPIInt                             mpiRank;


    /** \brief number of local packed variables of all processes. */
    types::IntVec1D                         nLclAllProcs;


    /** \brief offsets of packed variables of all processes. */
    types::IntVec1D                         offsetsAllProcs;


    /**
     * \brief create DMComposite.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode createDmPack();


    /**
     * \brief create ISLocalToGlobals for bodies.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode createMappings();


    /**
     * \brief create a string for printing information.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode createInfoString();
};

} // end of namespace body
} // end of namespace petibm
