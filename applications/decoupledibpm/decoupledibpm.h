/**
 * \file decoupledibpm.h
 * \brief Definition of the class \c DecoupledIBPMSolver.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 * \see decoupledibpm
 * \ingroup decoupledibpm
 */

#pragma once

#include <petibm/bodypack.h>

#include "../navierstokes/navierstokes.h"

/**
 * \class DecoupledIBPMSolver
 * \brief Immersed-boundary method proposed by Li et. al. (2016).
 * \see decoupledibpm, NavierStokesSolver
 * \ingroup decoupledibpm
 */
class DecoupledIBPMSolver : protected NavierStokesSolver
{
public:
    /** \brief Default constructor. */
    DecoupledIBPMSolver() = default;

    /** \brief Constructor; Initialize the decoupled IBPM solver.
     *
     * \param world [in] MPI communicator
     * \param node [in] YAML configuration settings
     */
    DecoupledIBPMSolver(const MPI_Comm &world, const YAML::Node &node);

    /** \brief Default destructor. */
    ~DecoupledIBPMSolver();

    /** \brief manually destroy data. */
    PetscErrorCode destroy();

    /** \brief Initialize the decoupled IBPM solver.
     *
     * \param world [in] MPI communicator
     * \param node [in] YAML configuration settings
     */
    PetscErrorCode init(const MPI_Comm &world, const YAML::Node &node);

    using NavierStokesSolver::ioInitialData;

    /** \brief Advance the solution by one time step. */
    PetscErrorCode advance();

    /** \brief Write solution and solver info to files. */
    PetscErrorCode write();

    using NavierStokesSolver::finished;

protected:
    /** \brief Pack of immersed bodies. */
    petibm::type::BodyPack bodies;

    /** \brief Linear solver for the Lagrangian forces. */
    petibm::type::LinSolver fSolver;

    /** \brief Spreading operator. */
    Mat H;

    /** \brief Regularization operator. */
    Mat E;

    /** \brief Left-hand side operator of the system for the forces. */
    Mat EBNH;

    /** \brief Projection operator for the forces. */
    Mat BNH;

    /** \brief Right-hand-side vector of the system for the forces. */
    Vec Eu;

    /** \brief Vector to hold the forces at time step n. */
    Vec f;

    /** \brief Force-increment vector. */
    Vec df;

    /** \brief Log stage for assembling the RHS of the forces system. */
    PetscLogStage stageRHSForces;

    /** \brief Log stage for solving the forces system. */
    PetscLogStage stageSolveForces;

    /** \brief Log stage for integrating the Lagrangian forces. */
    PetscLogStage stageIntegrateForces;

    /** \brief ASCII PetscViewer object to output the forces. */
    PetscViewer forcesViewer;

    /** \brief Assemble the RHS vector of the velocity system. */
    virtual PetscErrorCode assembleRHSVelocity();

    /** \brief Assemble the RHS vector of the Poisson system. */
    virtual PetscErrorCode assembleRHSPoisson();

    /** \brief Assemble the RHS vector of the system for the Lagrangian forces. */
    virtual PetscErrorCode assembleRHSForces();

    /** \brief Solve the system for the boundary forces. */
    virtual PetscErrorCode solveForces();

    /** \brief Project the velocity field, update pressure and forces. */
    virtual PetscErrorCode projectionStep();

    /** \brief Assemble additional operators. */
    PetscErrorCode createExtraOperators();

    /** \brief Create additional vectors. */
    PetscErrorCode createExtraVectors();

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

    /** \brief Write numbers of iterations and residuals of solvers to file. */
    PetscErrorCode writeLinSolversInfo();

    /** \brief Write the forces acting on the bodies into an ASCII file. */
    PetscErrorCode writeForcesASCII();

};  // DecoupledIBPMSolver
