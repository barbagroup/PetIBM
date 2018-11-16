#include <iomanip>

#include "rigidkinematics.h"

RigidKinematicsSolver::RigidKinematicsSolver(const MPI_Comm &world,
                                             const YAML::Node &node)
{
    init(world, node);
}  // RigidKinematicsSolver::RigidKinematicsSolver

RigidKinematicsSolver::~RigidKinematicsSolver()
{
    PetscErrorCode ierr;
    PetscBool finalized;

    PetscFunctionBeginUser;

    ierr = PetscFinalized(&finalized); CHKERRV(ierr);
    if (finalized) return;

    ierr = destroy(); CHKERRV(ierr);
}  // ~RigidKinematicsSolver::RigidKinematicsSolver

PetscErrorCode RigidKinematicsSolver::destroy()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = DecoupledIBPMSolver::destroy(); CHKERRQ(ierr);
    ierr = VecDestroy(&UB); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // RigidKinematicsSolver::destroy

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
}  // RigidKinematicsSolver::init

PetscErrorCode RigidKinematicsSolver::advance()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = moveBodies(t + dt); CHKERRQ(ierr);

    ierr = DecoupledIBPMSolver::advance(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // RigidKinematicsSolver::advance

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
}  // RigidKinematicsSolver::write

PetscErrorCode RigidKinematicsSolver::ioInitialData()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = DecoupledIBPMSolver::ioInitialData(); CHKERRQ(ierr);

    ierr = writeBodies(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // RigidKinematicsSolver

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
}  // RigidKinematicsSolver::moveBodies

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
}  // RigidKinematicsSolver::assembleRHSForces

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
}  // RigidKinematicsSolver::writeBodies
