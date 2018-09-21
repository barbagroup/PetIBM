/**
 * \file singlebody.h
 * \brief body::SingleBodyBase, type::SingleBody factory function.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 */

#pragma once

#include <memory>

#include <petscdm.h>
#include <petscsys.h>

#include <petibm/mesh.h>
#include <petibm/type.h>

namespace petibm
{
namespace body
{
// TODO: the way we store coords and meshIdx is different from CartesianMesh.
//       Will this cause confusion?
/**
 * \brief Base (abstract) class for a single body.
 * \see bodyModule, petibm::type::SingleBody, petibm::body::createSingleBody
 * \ingroup bodyModule
 *
 * There is currently only one implementation for this abstract class:
 * body::SingleBodyPoints.
 */
class SingleBodyBase
{
    friend class BodyPackBase;

public:
    /** \brief Number of dimensions. */
    PetscInt dim;

    /** \brief Name of the body. */
    std::string name;

    /** \brief Path of the file with the body coordinates. */
    std::string filePath;

    /** \brief Total number of Lagrangian points. */
    PetscInt nPts;

    /** \brief Coordinates of ALL Lagrangian points. */
    type::RealVec2D coords;

    /** \brief Local number of Lagrangian points. */
    PetscInt nLclPts;

    /**
     * \brief Index of the closest Eulerian mesh cell
     *        for each local Lagrangian point.
     */
    type::IntVec2D meshIdx;

    /** \brief Global index of the first local Lagrangian point. */
    PetscInt bgPt;

    /** \brief Global index of the last local Lagrangian point. */
    PetscInt edPt;

    /** \brief Parallel layout of the body (as a 1D DMDA object). */
    DM da;

    /** \brief String with information about the body. */
    std::string info;

    /**
     * \brief Constructor. Initialize the body.
     *
     * \param comm [in] MPI communicator.
     * \param dim [in] Number of dimensions.
     * \param name [in] Name of body.
     * \param filePath [in] Path of the file with data of the body.
     *
     * Data in `filePath` depend on the class used.
     * For example, with the class `body::SingleBodyPoints`, the `filePath`
     * should contain the number of Lagrangian points followed by the list of
     * coordinates (placed in columns, the first column being the x coordinate).
     * Note: `body::SingleBodyPoints` is currently the only sub-class
     * implemented.
     */
    SingleBodyBase(const MPI_Comm &comm, const PetscInt &dim,
                   const std::string &name, const std::string &filePath);

    /** \brief Default constructor. */
    SingleBodyBase() = default;

    /** \brief Default destructor. */
    virtual ~SingleBodyBase();

    /** \brief Manually destroy data. */
    virtual PetscErrorCode destroy();

    /**
     * \brief Print information about the body.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode printInfo() const;

    /**
     * \brief Get the process id owning a given Lagrangian point.
     *
     * \param i [in] Index of the Lagrangian point.
     * \param p [out] Process id.
     *
     * \return PetscErrorCode.
     *
     * Note: there is no need to provide which degree of freedom is requested
     * as all degrees for a Lagrangian point are on the same process.
     */
    virtual PetscErrorCode findProc(const PetscInt &i,
                                    PetscMPIInt &p) const = 0;

    /**
     * \brief Get the global index of a Lagrangian point in a DMDA object
     *        given its degree of freedom.
     *
     * \param i [in] Index of the Lagrangian point.
     * \param dof [in] Degree of freedom.
     * \param idx [out] Global index of the Lagrangian point.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode getGlobalIndex(const PetscInt &i,
                                          const PetscInt &dof,
                                          PetscInt &idx) const = 0;

    /**
     * \brief Get the global index of a Lagrangian point in a DMDA object
     *        given a MatStencil.
     *
     * \param s [in] MatStencil of the Lagrangian point.
     * \param idx [out] Global index of the Lagrangian point.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode getGlobalIndex(const MatStencil &s,
                                          PetscInt &idx) const = 0;

    /**
     * \brief Calculate the averaged force of this body.
     *
     * \param f [in] Vec of forces on each Lagrangian point of this body.
     * \param fAvg [out] return averaged force with length equal to dimension.
     *
     * \return PetscErrorCode.
     *
     * Note: fAvg should have the correct size.
     * This function does not check if fAvg has been allocated correctly.
     */
    virtual PetscErrorCode calculateAvgForces(const Vec &f,
                                              type::RealVec1D &fAvg) const = 0;

    /** \brief Get the index of closest Eulerian mesh cell
     *         for each local Lagrangian point.
     *
     * \param mesh [in] Structured Cartesian mesh.
     *
     * \return PetscErrorCode.
     */
    virtual PetscErrorCode updateMeshIdx(const type::Mesh &mesh) = 0;

protected:
    /** \brief MPI communicator. */
    MPI_Comm comm;

    /** \brief Total number of processes. */
    PetscMPIInt mpiSize;

    /** \brief Rank of the local process. */
    PetscMPIInt mpiRank;

    /** \brief Vector with the number of local unknowns on each process. */
    type::IntVec1D nLclAllProcs;

    /** \brief Offset on each process. */
    type::IntVec1D offsetsAllProcs;

};  // SingleBodyBase

}  // end of namespace body

namespace type
{
/**
 * \brief Definition of type::SingleBody.
 *
 * \see bodyModule, petibm::body::SingleBodyBase, petibm::body::createSingleBody
 * \ingroup bodyModule
 */
typedef std::shared_ptr<body::SingleBodyBase> SingleBody;
}  // end of namespace type

namespace body
{
/**
 * \brief Factory function to create a single body.
 *
 * \param comm [in] MPI communicator.
 * \param dim [in] Number of dimensions.
 * \param type [in] Type of file (currently only accept "point").
 * \param name [in] Name of the body.
 * \param filePath [in] Path of the file with the coordinates.
 * \param bodies [out] SingleBody data object.
 *
 * \return PetscErrorCode.
 *
 * \see bodyModule, petibm::type::SingleBody
 * \ingroup bodyModule
 */
PetscErrorCode createSingleBody(const MPI_Comm &comm, const PetscInt &dim,
                                const std::string &type,
                                const std::string &name,
                                const std::string &filePath,
                                type::SingleBody &body);

}  // end of namespace body

}  // end of namespace petibm
