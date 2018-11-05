#pragma once

#include "../decoupledibpm/decoupledibpm.h"

class RigidKinematics : protected DecoupledIBPMSolver
{
public:
    RigidKinematics() = default;
    RigidKinematics(const MPI_Comm &world, const YAML::Node &node);
    ~RigidKinematics();
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
        const PetscReal &time){PetscFunctionReturn(0);};
    PetscErrorCode assembleRHSForces();
    PetscErrorCode writeBodiesCoordinates();
    PetscErrorCode writeBodyCoordinates(const std::string filepath,
                                        const petibm::type::SingleBody body);

};  // RigidKinematics
