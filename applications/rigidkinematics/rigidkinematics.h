#pragma once

#include "../decoupledibpm/decoupledibpm.h"

class RigidKinematicsSolver : protected DecoupledIBPMSolver
{
public:
    RigidKinematicsSolver() = default;
    RigidKinematicsSolver(const MPI_Comm &world, const YAML::Node &node);
    ~RigidKinematicsSolver();
    PetscErrorCode destroy();
    PetscErrorCode init(const MPI_Comm &world, const YAML::Node &node);
    PetscErrorCode advance();
    PetscErrorCode write();
    PetscErrorCode ioInitialData();
    using DecoupledIBPMSolver::finished;

protected:
    Vec UB;
    PetscLogStage stageMoveIB;
    PetscErrorCode moveBodies(const PetscReal &ti);
    virtual PetscErrorCode setCoordinatesBodies(
        const PetscReal &ti){PetscFunctionReturn(0);};
    virtual PetscErrorCode setVelocityBodies(
        const PetscReal &ti){PetscFunctionReturn(0);};
    PetscErrorCode assembleRHSForces();
    PetscErrorCode writeBodies();

};  // RigidKinematicsSolver
