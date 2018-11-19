/**
 * \file rigidkinematics.cpp
 * \brief Implementation of the class \c RigidKinematicsSolver.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 * \see decoupledibpm
 * \ingroup decoupledibpm
 */

#include <iomanip>

#include "rigidkinematics.h"

RigidKinematicsSolver::RigidKinematicsSolver(const MPI_Comm &world,
                                             const YAML::Node &node)
{
    init(world, node);
}  // RigidKinematicsSolver

RigidKinematicsSolver::~RigidKinematicsSolver()
{
    PetscErrorCode ierr;
    PetscBool finalized;

    PetscFunctionBeginUser;

    ierr = PetscFinalized(&finalized); CHKERRV(ierr);
    if (finalized) return;

    ierr = destroy(); CHKERRV(ierr);
}  // ~RigidKinematicsSolver


// destroy data
PetscErrorCode RigidKinematicsSolver::destroy()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = DecoupledIBPMSolver::destroy(); CHKERRQ(ierr);
    ierr = VecDestroy(&UB); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // destroy

// initialize the decoupled IBPM solver for moving rigid bodies
PetscErrorCode RigidKinematicsSolver::init(const MPI_Comm &world,
                                           const YAML::Node &node)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = DecoupledIBPMSolver::init(world, node); CHKERRQ(ierr);

    ierr = PetscLogStagePush(stageInitialize); CHKERRQ(ierr);

    ierr = PetscLogStageRegister("moveIB", &stageMoveIB); CHKERRQ(ierr);

    ierr = VecDuplicate(f, &UB); CHKERRQ(ierr);
    ierr = VecSet(UB, 0.0); CHKERRQ(ierr);

    ierr = PetscLogStagePop(); CHKERRQ(ierr);  // end of stageInitialize

    PetscFunctionReturn(0);
}  // init

// advance the solution by one time step
PetscErrorCode RigidKinematicsSolver::advance()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    // note: use of `t + dt` because `t` is update in the Navier-Stokes method
    ierr = moveBodies(t + dt); CHKERRQ(ierr);

    ierr = DecoupledIBPMSolver::advance(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // advance

// write the solution, solver info, and body points to files
PetscErrorCode RigidKinematicsSolver::write()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = DecoupledIBPMSolver::write(); CHKERRQ(ierr);

    if (ite % nsave == 0)
    {
        ierr = PetscLogStagePush(stageWrite); CHKERRQ(ierr);

        ierr = writeBodies(); CHKERRQ(ierr);

        ierr = PetscLogStagePop(); CHKERRQ(ierr);  // end of stageWrite
    }

    PetscFunctionReturn(0);
}  // write

// read of write initial data
PetscErrorCode RigidKinematicsSolver::ioInitialData()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = DecoupledIBPMSolver::ioInitialData(); CHKERRQ(ierr);

    ierr = writeBodies(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // ioInitialData

// update Lagrangian points, boundary velocity, and operators
PetscErrorCode RigidKinematicsSolver::moveBodies(const PetscReal &ti)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = PetscLogStagePush(stageMoveIB); CHKERRQ(ierr);

    ierr = setCoordinatesBodies(ti); CHKERRQ(ierr);
    ierr = setVelocityBodies(ti); CHKERRQ(ierr);
    ierr = bodies->updateMeshIdx(mesh); CHKERRQ(ierr);
    if (E != PETSC_NULL) {ierr = MatDestroy(&E); CHKERRQ(ierr);}
    if (H != PETSC_NULL) {ierr = MatDestroy(&H); CHKERRQ(ierr);}
    if (BNH != PETSC_NULL) {ierr = MatDestroy(&BNH); CHKERRQ(ierr);}
    if (EBNH != PETSC_NULL) {ierr = MatDestroy(&EBNH); CHKERRQ(ierr);}
    ierr = createExtraOperators(); CHKERRQ(ierr);
    ierr = fSolver->setMatrix(EBNH); CHKERRQ(ierr);

    ierr = PetscLogStagePop(); CHKERRQ(ierr);  // end of stageMoveIB

    PetscFunctionReturn(0);
}  // moveBodies

// assemble the right-hand side of the system for the Lagrangian forces
PetscErrorCode RigidKinematicsSolver::assembleRHSForces()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = PetscLogStagePush(stageRHSForces); CHKERRQ(ierr);

    // rhsf = UB - E u^{**}
    ierr = MatMult(E, solution->UGlobal, rhsf); CHKERRQ(ierr);
    ierr = VecScale(rhsf, -1.0); CHKERRQ(ierr);
    ierr = VecAYPX(rhsf, 1.0, UB); CHKERRQ(ierr);

    ierr = PetscLogStagePop(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // assembleRHSForces

// write the coordinates of the Lagrangian points into files
PetscErrorCode RigidKinematicsSolver::writeBodies()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    std::string directory = config["output"].as<std::string>(".");
    std::stringstream ss;
    ss << std::setfill('0') << std::setw(7) << ite << ".3D";
    std::string filepath;
    for (PetscInt i = 0; i < bodies->nBodies; ++i)
    {
        petibm::type::SingleBody &body = bodies->bodies[i];
        filepath = directory + "/" + body->name + "_" + ss.str();
        ierr = body->writeBody(filepath); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}  // writeBodies
