/**
 * \file singlebodypoints.h
 * \brief Definition of body::SingleBodyPoints.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 */

#pragma once

#include <petibm/singlebody.h>

namespace petibm
{
namespace body
{
/**
 * \brief An implementation of body::SingleBodyBase that uses point data as
 * input.
 *
 * \see bodyModule, petibm::type::SingleBody, petibm::body::SingleBodyBase
 * \ingroup bodyModule
 *
 * This implementation uses an ASCII file of coordinates of Lagrangian points as
 * its input.
 *
 * Users should not initialize an instance of this class directly. They should
 * use petibm::body::createSingleBody. Actually, the design of PetIBM is to
 * handle multiple bodies, users should use petibm::type::BodyPack to handle
 * their bodies even if there is only one body.
 */
class SingleBodyPoints : public SingleBodyBase
{
public:
    /**
     * \brief Constructor. Initialize a single body.
     *
     * \param comm [in] MPI communicator.
     * \param dim [in] Number of dimensions.
     * \param name [in] Name of the body.
     * \param filePath [in] Path of the file with the coordinates.
     */
    SingleBodyPoints(const MPI_Comm &comm, const PetscInt &dim,
                     const std::string &name, const std::string &filePath);

    /** \copydoc SingleBodyBase::~SingleBodyBase */
    virtual ~SingleBodyPoints() = default;

    // implementation of SingleBodyBase::findProc
    virtual PetscErrorCode findProc(const PetscInt &i, PetscMPIInt &p) const;

    // implementation of SingleBodyBase::getGlobalIndex
    virtual PetscErrorCode getGlobalIndex(const PetscInt &i,
                                          const PetscInt &dof,
                                          PetscInt &idx) const;

    // implementation of SingleBodyBase::getGlobalIndex
    virtual PetscErrorCode getGlobalIndex(const MatStencil &s,
                                          PetscInt &idx) const;

    // implementation of SingleBodyBase::calculateAvgForces
    virtual PetscErrorCode calculateAvgForces(const Vec &f,
                                              type::RealVec1D &fAvg) const;

    // implementation of SingleBodyBase::updateMeshIdx
    virtual PetscErrorCode updateMeshIdx(const type::Mesh &mesh);

    virtual PetscErrorCode readBody(const std::string &filepath);

    virtual PetscErrorCode writeBody(const std::string &filepath);

protected:
    /**
     * \brief Initialize the body, reading coordinates from given file.
     *
     * \param comm [in] MPI communicator.
     * \param dim [in] Number of dimensions.
     * \param name [in] Name of the body.
     * \param filePath [in] Path of the file with the coordinates.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode init(const MPI_Comm &comm, const PetscInt &dim,
                        const std::string &name, const std::string &filePath);

    /**
     * \brief Create a parallel layout (1D DMDA object) of the body.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode createDMDA();

    /**
     * \brief Create a string used to print information about the body.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode createInfoString();

};  // SingleBodyPoints

}  // end of namespace body

}  // end of namespace petibm
