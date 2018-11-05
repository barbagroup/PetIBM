#include <iomanip>

#include "rigidkinematics.h"

RigidKinematics::RigidKinematics(const MPI_Comm &world,
                                 const YAML::Node &node)
{
    init(world, node);
}  // RigidKinematics

RigidKinematics::~RigidKinematics()
{
    PetscErrorCode ierr;
    PetscBool finalized;

    PetscFunctionBeginUser;

    ierr = PetscFinalized(&finalized); CHKERRV(ierr);
    if (finalized) return;

    ierr = destroy(); CHKERRV(ierr);
}  // ~RigidKinematics

PetscErrorCode RigidKinematics::destroy()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = DecoupledIBPMSolver::destroy(); CHKERRQ(ierr);
    ierr = VecDestroy(&UB); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // RigidKinematics::destroy

PetscErrorCode RigidKinematics::init(const MPI_Comm &world,
                                     const YAML::Node &node)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = DecoupledIBPMSolver::init(world, node); CHKERRQ(ierr);

    ierr = PetscLogStagePush(stageInitialize); CHKERRQ(ierr);

    ierr = PetscLogStageRegister("moveIB", &stageMoveIB); CHKERRQ(ierr);

    ierr = VecDuplicate(f, &UB); CHKERRQ(ierr);

    ierr = PetscLogStagePop(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // RigidKinematics::init

PetscErrorCode RigidKinematics::advance()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = PetscLogStagePush(stageMoveIB); CHKERRQ(ierr);
    ierr = moveIB(t + dt); CHKERRQ(ierr);
    ierr = createExtraOperators(); CHKERRQ(ierr);
    ierr = fSolver->setMatrix(EBNH); CHKERRQ(ierr);
    ierr = PetscLogStagePop(); CHKERRQ(ierr);

    ierr = DecoupledIBPMSolver::advance(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // RigidKinematics::advance

PetscErrorCode RigidKinematics::write()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = DecoupledIBPMSolver::write(); CHKERRQ(ierr);

    ierr = PetscLogStagePush(stageWrite); CHKERRQ(ierr);
    ierr = writeBodiesCoordinates(); CHKERRQ(ierr);
    ierr = PetscLogStagePop(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // write

PetscErrorCode RigidKinematics::assembleRHSForces()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = PetscLogStagePush(stageRHSForces); CHKERRQ(ierr);

    // rhsf = UB - E u^{**}
    ierr = MatMult(E, solution->UGlobal, rhsf); CHKERRQ(ierr);
    ierr = VecScale(rhsf, -1.0); CHKERRQ(ierr);
    ierr = VecAYPX(rhsf, 1.0, UB); CHKERRQ(ierr);

    // we choose to set the PETSc Vec object df with zeros
    ierr = VecSet(df, 0.0); CHKERRQ(ierr);

    ierr = PetscLogStagePop(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // RigidKinematics::assembleRHSForces

PetscErrorCode RigidKinematics::writeBodiesCoordinates()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    std::string directory = config["output"].as<std::string>(".");
    std::stringstream ss;
    ss << "/" << std::setfill('0') << std::setw(7) << ite << ".3D";
    std::string filepath;
    for (PetscInt i = 0; i < bodies->nBodies; ++i)
    {
        filepath = directory + "/body" + std::to_string(i) + ss.str();
        ierr = writeBodyCoordinates(
            filepath, (bodies->bodies)[0]); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}  // RigidKinematics::writeBodiesCoordinates

PetscErrorCode RigidKinematics::writeBodyCoordinates(
    const std::string filepath, const petibm::type::SingleBody body)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    PetscViewer viewer;
    ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr);
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
} // RigidKinematics::writeBodyCoordinates