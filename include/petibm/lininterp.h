/**
 * \file lininterp.h
 * \brief Prototype of the linear interpolation classes,
 *        definition of type::LinInterp, and factory function.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 */

#pragma once

#include <petscksp.h>
#include <petscsys.h>

#include <petibm/mesh.h>
#include <petibm/type.h>

namespace petibm
{
namespace misc
{
/**
 * \brief Abstract Base Class of a linear interpolation object.
 *
 * \see miscModule, petibm::type::LinInterp, petibm::misc::createLinInterp
 * \ingroup miscModule
 */
class LinInterpBase
{
public:
    /** \brief Default constructor. */
    LinInterpBase() = default;

    /** \brief Destructor. */
    ~LinInterpBase();

    /** \brief Manually destroy data in the object. */
    PetscErrorCode destroy();

    /** \brief Return the linearly interpolated value.
     *
     * \param da [in] Parallel layout of the field to interpolate
     * \param vec [in] Vector with the field data
     * \param val [out] Interpolated value
     * \return PetscErrorCode
     */
    PetscErrorCode getValue(const DM &da, const Vec &vec, PetscReal &val);

protected:
    /** \brief Coordinates of the point to interpolate. */
    type::RealVec1D target;

    /** \brief Coordinates of the front-bottom-left neighbor. */
    type::RealVec1D bl;

    /** \brief Coordinates of the back-top-right neighbor. */
    type::RealVec1D tr;

    /** \brief Gridline indices of the front-bottom-left neighbor. */
    type::IntVec1D idxDirs;

    /** \brief Container to hold the linear base for the target. */
    type::RealVec1D base_a;

    /** \brief Container to hold the elements of the interpolation matrix. */
    type::RealVec1D Op_a;

    /** \brief PETSc Vec object with neighbor values. */
    Vec sub;

    /** \brief PETSc Vec object with the linear base. */
    Vec base;

    /** \brief PETSc Vec object with interpolation coefficients. */
    Vec coeffs;

    /** \brief PETSc Mat object for the interpolation operator. */
    Mat Op;

    /** \brief PETSc KSP object used as a direct solver. */
    KSP ksp;

    /** \brief MPI communicator. */
    MPI_Comm comm;

    /** \brief Size of the MPI communicator. */
    PetscMPIInt commSize;

    /** \brief Rank of the local process in the communicator. */
    PetscMPIInt commRank;

    /** \brief Initialize the linear interpolation object.
     *
     * \param comm [in] MPI communicator
     * \param point [in] Coordinates of the target point to interpolate
     * \param mesh [in] Cartesian mesh object
     * \param field [in] Index of the field to monitor
     * \return PetscErrorCode
     */
    virtual PetscErrorCode init(const MPI_Comm &comm,
                                const type::RealVec1D &point,
                                const type::Mesh &mesh,
                                const type::Field &field) = 0;

    /** \brief Get the gridline indices of the front-bottom-left neighbor.
     *
     * \param mesh [in] Cartesian mesh object
     * \param field [in] Index of the field to monitor
     * \return PetscErrorCode
     */
    PetscErrorCode getBLGridlineIndices(const type::Mesh &mesh,
                                        const type::Field &field);

    /** \brief Get the coordinates of the neighbors.
     *
     * \param mesh [in] Cartesian mesh object
     * \param field [in] Index of the field to monitor
     * \return PetscErrorCode
     */
    PetscErrorCode getBoxCoords(const type::Mesh &mesh,
                                const type::Field &field);

    /** \brief Set up the linear base. */
    virtual PetscErrorCode setUpBase() = 0;

    /** \brief Set up the linear interpolation operator. */
    virtual PetscErrorCode setUpOp() = 0;

    /** \brief Set up the direct solver
     *         to find the interpolation coefficients.
     */
    PetscErrorCode setUpKSP();

    /** \brief Set up the vector holding the neighbor values.
     *
     * \param da [in] PETSc DMDA object of the field to monitor
     * \param vec [in] PETSc Vec object with full-domain data
     * \return PetscErrorCode
     */
    virtual PetscErrorCode setSubVec(const DM &da, const Vec &vec) = 0;

};  // LinInterpBase

/** \brief Class to perform a trilinear interpolation at a point in a 3D domain.
 *
 * \see miscModule, petibm::type::LinInterp, petibm::misc::createLinInterp
 * \ingroup miscModule
 */
class TriLinInterp: public LinInterpBase
{
public:
    /** \brief Constructor. Initialize the linear interpolation object.
     *
     * \param comm [in] MPI communicator
     * \param point [in] Coordinates of the target point to interpolate
     * \param mesh [in] Cartesian mesh object
     * \param field [in] Index of the field to monitor
     * \return PetscErrorCode
     */
    TriLinInterp(const MPI_Comm &comm,
                 const type::RealVec1D &point,
                 const type::Mesh &mesh,
                 const type::Field &field);
    
    /** \brief Default destructor. */
    ~TriLinInterp() = default;

protected:
    /** copydoc LinInterpBase::LinInterpBase() */
    PetscErrorCode init(const MPI_Comm &comm,
                        const type::RealVec1D &point,
                        const type::Mesh &mesh,
                        const type::Field &field);

    /** copydoc LinInterpBase::setUpBase() */
    PetscErrorCode setUpBase();

    /** copydoc LinInterpBase::setUpOp() */
    PetscErrorCode setUpOp();

    /** copydoc LinInterpBase::setSubVec() */
    PetscErrorCode setSubVec(const DM &da, const Vec &vec);

};  // TriLinInterp

/** \brief Class to perform a bilinear interpolation at a point in a 3D domain.
 *
 * \see miscModule, petibm::type::LinInterp, petibm::misc::createLinInterp
 * \ingroup miscModule
 */
class BiLinInterp: public LinInterpBase
{
public:
    /** \brief Constructor. Initialize the linear interpolation object.
     *
     * \param comm [in] MPI communicator
     * \param point [in] Coordinates of the target point to interpolate
     * \param mesh [in] Cartesian mesh object
     * \param field [in] Index of the field to monitor
     * \return PetscErrorCode
     */
    BiLinInterp(const MPI_Comm &comm,
                const type::RealVec1D &point,
                const type::Mesh &mesh,
                const type::Field &field);
    
    /** \brief Default destructor. */
    ~BiLinInterp() = default;

protected:
    /** copydoc LinInterpBase::LinInterpBase() */
    PetscErrorCode init(const MPI_Comm &comm,
                        const type::RealVec1D &point,
                        const type::Mesh &mesh,
                        const type::Field &field);

    /** copydoc LinInterpBase::setUpBase() */
    PetscErrorCode setUpBase();

    /** copydoc LinInterpBase::setUpOp() */
    PetscErrorCode setUpOp();

    /** copydoc LinInterpBase::setSubVec() */
    PetscErrorCode setSubVec(const DM &da, const Vec &vec);

};  // BiLinInterp

}  // end of namespace misc

namespace type
{
/** \brief Type definition of LinInterp.
 *
 * Please use petibm::misc::createLinInterp to create a LinInterp object.
 *
 * * \see miscModule, petibm::misc::createLinInterp, petibm::misc::LinInterpBase
 * \ingroup miscModule
 */
typedef std::shared_ptr<misc::LinInterpBase> LinInterp;

}  // end of namespace type

namespace misc
{
/**
 * \brief Factory function to create a linear interpolation object.
 *
 * \param comm [in] MPI communicator
 * \param point [in] Coordinates of the point to interpolate
 * \param mesh [in] Cartesian mesh object
 * \param field [in] Index of the field to monitor
 * \param interp [out] Linear interpolation object
 * \return PetscErrorCode
 *
 * \see miscModule, petibm::type::LinInterp
 * \ingroup miscModule
 */
PetscErrorCode createLinInterp(const MPI_Comm &comm,
                               const type::RealVec1D &point,
                               const type::Mesh &mesh,
                               const type::Field &field,
                               type::LinInterp &interp);

}  // end of namespace misc

}  // end of namespace petibm
