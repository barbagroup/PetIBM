/**
 * \file oscillatingcylinder.cpp
 * \brief Implementation of the class \c OscillatingCylinderSolver.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 */

#include "oscillatingcylinder.h"

#include <petibm/io.h>

OscillatingCylinderSolver::OscillatingCylinderSolver(const MPI_Comm &world, const YAML::Node & node)
{
    init(world, node);
}  // OscillatingCylinderSolver

OscillatingCylinderSolver::~OscillatingCylinderSolver()
{
    PetscErrorCode ierr;
    PetscBool finalized;

    PetscFunctionBeginUser;

    ierr = PetscFinalized(&finalized); CHKERRV(ierr);
    if (finalized) return;

    ierr = destroy(); CHKERRV(ierr);
}  // OscillatingCylinderSolver::~OscillatingCylinderSolver

PetscErrorCode OscillatingCylinderSolver::init(const MPI_Comm &world,
                                               const YAML::Node &node)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    // initialize the decoupled IBPM solver for rigid body motion
    ierr = RigidKinematicsSolver::init(world, node); CHKERRQ(ierr);

    ierr = PetscLogStagePush(stageInitialize); CHKERRQ(ierr);

    // parse the configuration file to get kinematics parameters
    const YAML::Node &config_kin = node["bodies"][0]["kinematics"];
    // frequency of the oscillation
    f = config_kin["f"].as<PetscReal>(0.0);
    // diameter of the cylinder
    PetscReal Dc = config_kin["D"].as<PetscReal>(1.0);
    // Keulegan-Carpenter number
    PetscReal KC = config_kin["KC"].as<PetscReal>(0.0);
    // amplitude of the oscillation
    Am = Dc * KC / (2.0 * PETSC_PI);
    // maximum translational velocity
    Um = 2.0 * PETSC_PI * f * Am;
    // Initial location of the cylinder center
    Xc0 = config_kin["center"][0].as<PetscReal>(0.0);
    Yc0 = config_kin["center"][1].as<PetscReal>(0.0);

    ierr = PetscLogStagePop(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // init

PetscErrorCode OscillatingCylinderSolver::setCoordinatesBodies(
    const PetscReal &ti)
{
    PetscReal Xd,  // displacement in the x-direction
              Yd;  // displacement in the y-direction
    // make references for code readability
    petibm::type::SingleBody &body = bodies->bodies[0];
    petibm::type::RealVec2D &coords = body->coords;
    petibm::type::RealVec2D &coords0 = body->coords0;

    PetscFunctionBeginUser;

    // compute the displacements in the x and y directions
    Xd = -Am * PetscSinReal(2*PETSC_PI * f * ti);
    Yd = 0.0;

    // update the position of the body points
    for (PetscInt k = 0; k < body->nPts; k++)
    {
        coords[k][0] = coords0[k][0] + Xd;
        coords[k][1] = coords0[k][1] + Yd;
    }

    PetscFunctionReturn(0);
} // setCoordinatesBodies

PetscErrorCode OscillatingCylinderSolver::setVelocityBodies(
    const PetscReal &ti)
{
    PetscErrorCode ierr;
    PetscReal Ux;  // translation velocity in x-direction
    PetscReal **UB_arr;
    petibm::type::SingleBody &body = bodies->bodies[0];

    PetscFunctionBeginUser;

    // compute the translational velocity at current time
    Ux = - Um * PetscCosReal(2 * PETSC_PI * f * ti);
    // update the body velocity array
    ierr = DMDAVecGetArrayDOF(body->da, UB, &UB_arr); CHKERRQ(ierr);
    for (PetscInt k = body->bgPt; k < body->edPt; k++)
    {
        UB_arr[k][0] = Ux;
        UB_arr[k][1] = 0.0;
    }
    ierr = DMDAVecRestoreArrayDOF(body->da, UB, &UB_arr); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // setVelocityBodies
