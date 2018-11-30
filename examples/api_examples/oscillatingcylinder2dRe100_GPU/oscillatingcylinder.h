/**
 * \file oscillatingcylinder.h
 * \brief Definition of the class \c OscillatingCylinderSolver.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 */

#pragma once

#include <petibm/rigidkinematics/rigidkinematics.h>

/**
 * \class OscillatingCylinderSolver
 * \brief Solve the two-dimensional flow around an inline-oscillating cylinder
 *        with a provided rigid-body kinematics and using the decoupled IBPM.
 * \see decoupledibpm
 */
class OscillatingCylinderSolver : protected RigidKinematicsSolver
{
public:
    /** \brief Default constructor. */
    OscillatingCylinderSolver() = default;

    /** \brief Constructor; Initialize the decoupled IBPM solver with moving cylinder.
     *
     * \param world [in] MPI communicator
     * \param node [in] YAML configuration settings
     */
    OscillatingCylinderSolver(const MPI_Comm &world, const YAML::Node &node);

    /** \brief Default destructor. */
    ~OscillatingCylinderSolver();

    /** \brief Initialize the decoupled IBPM solver with moving rigid bodies.
     *
     * \param world [in] MPI communicator
     * \param node [in] YAML configuration settings
     */
    PetscErrorCode init(const MPI_Comm &world, const YAML::Node &node);

    /** \brief Manually destroy data. */
    using RigidKinematicsSolver::destroy;

    /** \brief Advance the solution by one time step. */
    using RigidKinematicsSolver::advance;

    /** \brief Write solution, solver info, and body points to files. */
    using RigidKinematicsSolver::write;

    /** \brief Read or write initial data. */
    using RigidKinematicsSolver::ioInitialData;

    /** \brief Check if the simulation is finished. */
    using RigidKinematicsSolver::finished;

protected:
    /** \brief Frequency of oscillation. */
    PetscReal f;

    /** \brief Amplitude of oscillation. */
    PetscReal Am;

    /** \brief Maximum translational velocity of the cylinder. */
    PetscReal Um;

    /** \brief Initial x-location of the cylinder center. */
    PetscReal Xc0;

    /** \brief Initial y-location of the cylinder center. */
    PetscReal Yc0;

    /** \brief Update the position of the Lagrangian points.
     *
     * \param ti [in] Time
     */
    PetscErrorCode setCoordinatesBodies(const PetscReal &ti);

    /** \brief Update the velocity of the Lagrangian points.
     *
     * \param ti [in] Time
     */
    PetscErrorCode setVelocityBodies(const PetscReal &ti);

};  // OscillatingCylinderSolver
