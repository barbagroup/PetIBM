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

    ierr = PetscLogStagePop(); CHKERRQ(ierr);  // end of stageInitialize

    PetscFunctionReturn(0);
}  // RigidKinematicsSolver::init

PetscErrorCode RigidKinematicsSolver::advance()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = PetscLogStagePush(stageMoveIB); CHKERRQ(ierr);

    ierr = moveIB(t + dt); CHKERRQ(ierr);
    ierr = bodies->updateMeshIdx(mesh); CHKERRQ(ierr);
    if (E != PETSC_NULL) {ierr = MatDestroy(&E); CHKERRQ(ierr);}
    if (H != PETSC_NULL) {ierr = MatDestroy(&H); CHKERRQ(ierr);}
    if (BNH != PETSC_NULL) {ierr = MatDestroy(&BNH); CHKERRQ(ierr);}
    if (EBNH != PETSC_NULL) {ierr = MatDestroy(&EBNH); CHKERRQ(ierr);}
    ierr = createExtraOperators(); CHKERRQ(ierr);
    ierr = fSolver->setMatrix(EBNH); CHKERRQ(ierr);

    ierr = PetscLogStagePop(); CHKERRQ(ierr);  // end of stageMoveIB

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

        ierr = writeBodiesPoint3D(); CHKERRQ(ierr);

        ierr = PetscLogStagePop(); CHKERRQ(ierr);  // end of stageWrite
    }

    PetscFunctionReturn(0);
}  // RigidKinematicsSolver::write

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

PetscErrorCode RigidKinematicsSolver::writeBodiesPoint3D()
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
        ierr = writeBodyPoint3D(filepath, body); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}  // RigidKinematicsSolver::writeBodiesPoint3D

PetscErrorCode RigidKinematicsSolver::writeBodyPoint3D(
    const std::string filepath, const petibm::type::SingleBody body)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    PetscViewer viewer;
    ierr = PetscViewerCreate(comm, &viewer); CHKERRQ(ierr);
    ierr = PetscViewerSetType(viewer, PETSCVIEWERASCII); CHKERRQ(ierr);
    ierr = PetscViewerFileSetMode(viewer, FILE_MODE_WRITE); CHKERRQ(ierr);
    ierr = PetscViewerFileSetName(viewer, filepath.c_str()); CHKERRQ(ierr);
    for (PetscInt k=0; k<body->nPts; k++)
    {
        ierr = PetscViewerASCIIPrintf(viewer, "%10.8e\t%10.8e\t%10.8e\n",
                                      body->coords[k][0],
                                      body->coords[k][1],
                                      body->coords[k][2]); CHKERRQ(ierr);
    }
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // RigidKinematicsSolver::writeBodyPoint3D
