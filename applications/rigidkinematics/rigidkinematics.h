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
    using DecoupledIBPMSolver::ioInitialData;
    using DecoupledIBPMSolver::finished;

protected:
    Vec UB;
    PetscLogStage stageMoveIB;
    virtual PetscErrorCode moveIB(
        const PetscReal &ti){PetscFunctionReturn(0);};
    PetscErrorCode assembleRHSForces();
    PetscErrorCode writeBodiesPoint3D();
    PetscErrorCode writeBodyPoint3D(const std::string filepath,
                                    const petibm::type::SingleBody body);

};  // RigidKinematicsSolver
