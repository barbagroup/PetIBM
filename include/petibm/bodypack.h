/**
 * \file bodypack.h
 * \brief body::BodyPackBase, type::BodyPack, and factory function.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include <petscdm.h>
#include <petscsys.h>
#include <yaml-cpp/yaml.h>

#include <petibm/singlebody.h>

/**
 * \defgroup bodyModule Immersed-boundary bodies
 * \brief Objects related to immersed-boundary bodies.
 *
 * This module contains objects related to immersed-boundary bodies and factory
 * functions.
 *
 * A petibm::type::SingleBody is an object holding information of a single body,
 * including coordinates of Lagrangian points, and the indices of relevant
 * background Eulerian mesh points, etc. It also has member functions regarding
 * to a single body. For example, to update the indices of relevant Eulerian
 * mesh points, I/O, calculating forces, etc. A petibm::type::SingleBody is
 * distributed to all MPI processes regardless the decomposition of the
 * background Eulerian mesh. It also handles the indexing of Lagrangian points
 * in a parallel PETSc Vec.
 *
 * A petibm::type::BodyPack is just a collection of all petibm::type::SingleBody
 * in the domain for the purpose of easy coding. A member function in a
 * petibm::type::BodyPack is generally just calling the corresponding function
 * inside each single body. petibm::type::BodyPack also works even there is only
 * one single body. It also handles the indexing of Lagrangian points in a
 * packed parallel PETSc Vec.
 *
 * Users should use petibm::body::createBodyPack to create an instance. Although
 * it's not encouraged, users can also use petibm::body::createSingleBody if
 * there is only one body in their application code.
 *
 * \see petibm::type::BodyPack, petibm::type::SingleBody,
 * petibm::body::createBodyPack \ingroup petibm
 */

namespace petibm
{
/**
 * \brief A namespace for objects and functions related to immersed-boundary
 * bodies. \see petibm::type::BodyPack, petibm::type::SingleBody,
 * petibm::body::createBodyPack \ingroup bodyModule
 */
namespace body
{
/**
 * \brief Base (abstract) class for a pack of bodies.
 * \see bodyModule, petibm::type::BodyPack, petibm::type::SingleBody,
 * petibm::body::createBodyPack \ingroup bodyModule
 *
 * This class is designed to be an abstract class though, it actually has full
 * implementations of all members because there is currently no necessary to do
 * otherwise. It's just simply a collection of bodies, though it also deals with
 * indexing of all Lagrangian points in a globally packed Vec.
 */
class BodyPackBase
{
public:
    /** \brief Dimension. */
    PetscInt dim;

    /** \brief Number of bodies in this pack. */
    PetscInt nBodies;

    /** \brief Total number of Lagrangian points. */
    PetscInt nPts;

    /** \brief Total number of local Lagrangian points. */
    PetscInt nLclPts;

    /** \brief A vector of SingleBody instances. */
    std::vector<type::SingleBody> bodies;

    /** \brief A DMComposite of DMs of all `SingleBody`s. */
    DM dmPack;

    /** \brief A string for printing information. */
    std::string info;

    /** \brief Default constructor. */
    BodyPackBase() = default;

    /**
     * \brief Constructor of using a type::Mesh and a YAML node.
     *
     * \param mesh [in] Strctured Cartesian mesh object.
     * \param node [in] YAML configurations.
     */
    BodyPackBase(const type::Mesh &mesh, const YAML::Node &node);

    /** \brief Default destructor. */
    virtual ~BodyPackBase();

    /** \brief Manually destroy data.  */
    virtual PetscErrorCode destroy();

    /** \brief Print information.  */
    PetscErrorCode printInfo() const;

    /**
     * \brief Find which process owns the target Lagrangian point of target
     * body.
     *
     * \param bIdx [in] Index of target body.
     * \param ptIdx [in] Index of target point.
     * \param proc [out] Process index.
     */
    PetscErrorCode findProc(const PetscInt &bIdx, const PetscInt &ptIdx,
                            PetscMPIInt &proc) const;

    /**
     * \brief Find unpacked global index of a DoF of Lagrangian point of a body.
     *
     * \param bIdx [in] Index of target body.
     * \param ptIdx [in] Index of target point.
     * \param dof [in] Index of target DoF.
     * \param idx [out] Unpacked global index.
     */
    PetscErrorCode getGlobalIndex(const PetscInt &bIdx, const PetscInt &ptIdx,
                                  const PetscInt &dof, PetscInt &idx) const;

    /**
     * \brief Find unpacked global index of a DoF of Lagrangian point of a body.
     *
     * \param bIdx [in] Index of target body.
     * \param s [in] PETSc MatStencil object of target point.
     * \param idx [out] Unpacked global index.
     */
    PetscErrorCode getGlobalIndex(const PetscInt &bIdx, const MatStencil &s,
                                  PetscInt &idx) const;

    /**
     * \brief Find packed global index of a DoF of Lagrangian point of a body.
     *
     * \param bIdx [in] Index of target body.
     * \param ptIdx [in] Index of target point.
     * \param dof [in] Index of target DoF.
     * \param idx [out] Packed global index.
     */
    PetscErrorCode getPackedGlobalIndex(const PetscInt &bIdx,
                                        const PetscInt &ptIdx,
                                        const PetscInt &dof,
                                        PetscInt &idx) const;

    /**
     * \brief Find packed global index of a DoF of Lagrangian point of a body.
     *
     * \param bIdx [in] Index of target body.
     * \param s [in] PETSc MatStencil object of target point.
     * \param idx [in] Packed global index.
     */
    PetscErrorCode getPackedGlobalIndex(const PetscInt &bIdx,
                                        const MatStencil &s,
                                        PetscInt &idx) const;

    /**
     * \brief Calculate the averaged force of each body.
     *
     * \param f [in] Packed force Vec of Lagrangian points.
     * \param fAvg [out] Averaged force for each body.
     *
     * Note: if `fAvg` doesn't have the correct size, the function will resize
     * it.
     */
    PetscErrorCode calculateAvgForces(const Vec &f, type::RealVec2D &fAvg);

    /** \brief Update the indices of background Pressure cells for all body. */
    PetscErrorCode updateMeshIdx();

protected:
    /**
     * \brief Constructor of using a type::Mesh and a YAML node.
     *
     * \param mesh [in] Structured Cartesian mesh object.
     * \param node [in] YAML configurations.
     */
    PetscErrorCode init(const type::Mesh &mesh, const YAML::Node &node);

    /** \brief Reference to background mesh. */
    type::Mesh mesh;

    /** \brief Reference to the MPI communicator. */
    MPI_Comm comm;

    /** \brief The total number of processes. */
    PetscMPIInt mpiSize;

    /** \brief The rank of this process. */
    PetscMPIInt mpiRank;

    /** \brief Number of local packed variables of all processes. */
    type::IntVec1D nLclAllProcs;

    /** \brief Offsets of packed variables of all processes. */
    type::IntVec1D offsetsAllProcs;

    /** \brief Create DMComposite.  */
    PetscErrorCode createDmPack();

    /** \brief Create a string for printing information.  */
    PetscErrorCode createInfoString();

};  // BodyPackBase

}  // end of namespace body

namespace type
{
/**
 * \brief Definition of type::BodyPack.
 * \see bodyModule, petibm::type::SingleBody, petibm::body::createBodyPack
 * \ingroup bodyModule
 */
typedef std::shared_ptr<body::BodyPackBase> BodyPack;

}  // end of namespace type

namespace body
{
/**
 * \brief Factory function for creating type::BodyPack
 *
 * \param mesh [in] Structured Cartesian mesh object.
 * \param node [in] YAML configurations.
 * \param bodies [out] Data object for the body pack.
 *
 * \see bodyModule, petibm::type::BodyPack, petibm::type::SingleBody
 * \ingroup bodyModule
 */
PetscErrorCode createBodyPack(const type::Mesh &mesh, const YAML::Node &node,
                              type::BodyPack &bodies);

}  // end of namespace body

}  // end of namespace petibm
