/**
 * \file rigidkinematics.h
 * \brief Definition of the class \c RigidKinematicsSolver.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 * \see decoupledibpm
 * \ingroup decoupledibpm
 */

#pragma once

#include "../decoupledibpm/decoupledibpm.h"

/**
 * \class RigidKinematicsSolver
 * \brief Helper class to solve cases with moving rigid bodies using the decoupled IBPM.
 * \see decoupledibpm
 * \ingroup decoupledibpm
 */
class RigidKinematicsSolver : protected DecoupledIBPMSolver
{
public:
    /** \brief Default constructor. */
    RigidKinematicsSolver() = default;

    /** \brief Constructor; Initialize the decoupled IBPM solver with moving bodies.
     *
     * \param world [in] MPI communicator
     * \param node [in] YAML configuration settings
     */
    RigidKinematicsSolver(const MPI_Comm &world, const YAML::Node &node);

    /** \brief Default destructor. */
    ~RigidKinematicsSolver();

    /** \brief Manually destroy data. */
    PetscErrorCode destroy();

    /** \brief Initialize the decoupled IBPM solver with moving rigid bodies.
     *
     * \param world [in] MPI communicator
     * \param node [in] YAML configuration settings
     */
    PetscErrorCode init(const MPI_Comm &world, const YAML::Node &node);

    /** \brief Advance the solution by one time step. */
    PetscErrorCode advance();

    /** \brief Write solution, solver info, and body points to files. */
    PetscErrorCode write();

    /** \brief Read or write initial data. */
    PetscErrorCode ioInitialData();

    using DecoupledIBPMSolver::finished;

protected:
    /** \brief Prescribed boundary velocity vector. */
    Vec UB;

    /** \brief Log stage for moving the bodies. */
    PetscLogStage stageMoveIB;

    /** \brief Move the bodies and update the solver state.
     *
     * Update the position of the Lagrangian points.
     * Set the velocity of the Lagrangian points.
     * Assemble the operators that depends on the location of the Lagrangian points.
     *
     * \params ti [in] Time
     */
    virtual PetscErrorCode moveBodies(const PetscReal &ti);

    /** \brief Update the position of the Lagrangian points.
     *
     * The present implementation does nothing.
     * The method should be implemented by the user in a child class.
     *
     * \params ti [in] Time
     */
    virtual PetscErrorCode setCoordinatesBodies(
        const PetscReal &ti){PetscFunctionReturn(0);};
    
    /** \brief Update the velocity of the Lagrangian points.
     *
     * The present implementation does nothing.
     * The method should be implemented by the user in a child class.
     *
     * \params ti [in] Time
     */
    virtual PetscErrorCode setVelocityBodies(
        const PetscReal &ti){PetscFunctionReturn(0);};
    
    /** \brief Assemble the right-hand side of the system for the forces. */
    virtual PetscErrorCode assembleRHSForces();

    /** \brief Write the coordinates of the Lagrangian points for each body. */
    virtual PetscErrorCode writeBodies();

};  // RigidKinematicsSolver
