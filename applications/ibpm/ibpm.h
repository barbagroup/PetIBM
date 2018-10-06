/**
 * \file ibpm.h
 * \brief Definition of the class \c IBPMSolver.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 * \see ibpm
 * \ingroup ibpm
 */

#pragma once

#include <petibm/bodypack.h>

#include "../navierstokes/navierstokes.h"

/**
 * \class IBPMSolver
 * \brief Immersed-boundary method proposed by Taira and Colonius (2007).
 * \see ibpm, NavierStokesSolver
 * \ingroup ibpm
 */
class IBPMSolver : protected NavierStokesSolver
{
public:
    /** \brief Default constructor. */
    IBPMSolver() = default;

    /**
     * \brief Constructor. Initialize the IBPM solver.
     *
     * \param world [in] MPI communicator
     * \param node [in] YAML configuration settings
     */
    IBPMSolver(const MPI_Comm &world, const YAML::Node &node);

    /** \brief Default destructor. */
    ~IBPMSolver();

    /** \brief Manually destroy data. */
    PetscErrorCode destroy();

    /** \brief Initialize the IBPM solver.
     *
     * \param world [in] MPI communicator
     * \param node [in] YAML configuration settings
     * \return PetscErrorCode
     */
    PetscErrorCode init(const MPI_Comm &world, const YAML::Node &node);

    using NavierStokesSolver::ioInitialData;

    using NavierStokesSolver::advance;

    /** \brief Write solution, forces, and solvers info to files. */
    PetscErrorCode write();

    using NavierStokesSolver::finished;

protected:
    /** \brief Pack of immersed bodies. */
    petibm::type::BodyPack bodies;

    /** \brief PETSc Vec object with pressure field and Lagrangian forces. */
    Vec phi;

    /** \brief Global index sets for pressure field and Lagrangian forces. */
    IS isDE[2];

    /** \brief Log stage for integrating the Lagrangian forces. */
    PetscLogStage stageIntegrateForces;

    /** \brief ASCII PetscViewer object to output the forces. */
    PetscViewer forcesViewer;

    /** \brief Assemble the RHS vector of the Poisson system. */
    virtual PetscErrorCode assembleRHSPoisson();

    /** \brief Create operators. */
    virtual PetscErrorCode createOperators();

    /** \brief Create vectors. */
    virtual PetscErrorCode createVectors();

    /** \brief Set Poisson nullspace or pin pressure at a reference point. */
    virtual PetscErrorCode setNullSpace();

    /** \brief Write the solution fields into a HDF5 file.
     *
     * \param filePath [in] Path of the file to write in
     * \return PetscErrorCode
     */
    PetscErrorCode writeSolutionHDF5(const std::string &filePath);

    /** \brief Write data required to restart a simulation into a HDF5 file.
     *
     * \param filePath [in] Path of the file to write in
     * \return PetscErrorCode
     */
    PetscErrorCode writeRestartDataHDF5(const std::string &filePath);

    /** \brief Read data required to restart a simulation from a HDF5 file.
     *
     * \param filePath [in] Path of the file to read from
     * \return PetscErrorCode
     */
    PetscErrorCode readRestartDataHDF5(const std::string &filePath);

    /** \brief Write the forces acting on the bodies into an ASCII file. */
    PetscErrorCode writeForcesASCII();

};  // IBPMSolver
