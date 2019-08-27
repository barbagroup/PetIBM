/**
 * \file decoupledibpm.cpp
 * \brief Implementation of the class \c DecoupledIBPMSolver.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 * \see decoupledibpm
 * \ingroup decoupledibpm
 */

#include <petscviewerhdf5.h>

#include <petibm/delta.h>

#include "decoupledibpm.h"

DecoupledIBPMSolver::DecoupledIBPMSolver(const MPI_Comm &world,
                                         const YAML::Node &node)
{
    init(world, node);
}  // DecoupledIBPMSolver

DecoupledIBPMSolver::~DecoupledIBPMSolver()
{
    PetscErrorCode ierr;
    PetscBool finalized;

    PetscFunctionBeginUser;

    ierr = PetscFinalized(&finalized); CHKERRV(ierr);
    if (finalized) return;

    ierr = destroy(); CHKERRV(ierr);
}  // ~DecoupledIBPMSolver

// destroy
PetscErrorCode DecoupledIBPMSolver::destroy()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    fSolver.reset();
    bodies.reset();
    ierr = VecDestroy(&df); CHKERRQ(ierr);
    ierr = VecDestroy(&f); CHKERRQ(ierr);
    ierr = VecDestroy(&rhsf); CHKERRQ(ierr);
    ierr = MatDestroy(&H); CHKERRQ(ierr);
    ierr = MatDestroy(&E); CHKERRQ(ierr);
    ierr = MatDestroy(&EBNH); CHKERRQ(ierr);
    ierr = MatDestroy(&BNH); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&forcesViewer); CHKERRQ(ierr);
    ierr = NavierStokesSolver::destroy(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // destroy

PetscErrorCode DecoupledIBPMSolver::init(const MPI_Comm &world,
                                         const YAML::Node &node)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = NavierStokesSolver::init(world, node); CHKERRQ(ierr);

    ierr = PetscLogStagePush(stageInitialize); CHKERRQ(ierr);

    // create a pack of immersed bodies
    ierr = petibm::body::createBodyPack(
        comm, mesh->dim, config, bodies); CHKERRQ(ierr);
    ierr = bodies->updateMeshIdx(mesh); CHKERRQ(ierr);

    // create the linear solver object for the Lagrangian forces
    ierr = petibm::linsolver::createLinSolver(
        "forces", config, fSolver); CHKERRQ(ierr);

    // create additional operators required for the decoupled IBPM
    ierr = createExtraOperators(); CHKERRQ(ierr);

    // create additional vectors required for the decoupled IBPM
    ierr = createExtraVectors(); CHKERRQ(ierr);

    // set coefficient matrix to the linear solver for the forces
    ierr = fSolver->setMatrix(EBNH); CHKERRQ(ierr);

    // create an ASCII PetscViewer to output the body forces
    ierr = createPetscViewerASCII(
        config["directory"].as<std::string>() +
        "/forces-" + std::to_string(ite) + ".txt",
        FILE_MODE_WRITE, forcesViewer); CHKERRQ(ierr);

    // register additional logging stages
    ierr = PetscLogStageRegister("rhsForces", &stageRHSForces); CHKERRQ(ierr);
    ierr = PetscLogStageRegister(
        "solveForces", &stageSolveForces); CHKERRQ(ierr);
    ierr = PetscLogStageRegister(
        "integrateForces", &stageIntegrateForces); CHKERRQ(ierr);

    ierr = PetscLogStagePop(); CHKERRQ(ierr);  // end of stageInitialize

    PetscFunctionReturn(0);
}  // init

// advance the decoupled IBPM solver by one time-step
PetscErrorCode DecoupledIBPMSolver::advance()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    
    t += dt;
    ite++;
    
    ierr = assembleRHSVelocity(); CHKERRQ(ierr);
    ierr = solveVelocity(); CHKERRQ(ierr);
    
    ierr = assembleRHSForces(); CHKERRQ(ierr);
    ierr = solveForces(); CHKERRQ(ierr);
    ierr = applyNoSlip(); CHKERRQ(ierr);

    ierr = assembleRHSPoisson(); CHKERRQ(ierr);
    ierr = solvePoisson(); CHKERRQ(ierr);
    ierr = applyDivergenceFreeVelocity(); CHKERRQ(ierr);

    ierr = updatePressure(); CHKERRQ(ierr);
    ierr = updateForces(); CHKERRQ(ierr);

    ierr = bc->updateGhostValues(solution); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // advance

// write solution fields, linear solvers info, and body forces to files
PetscErrorCode DecoupledIBPMSolver::write()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = NavierStokesSolver::write(); CHKERRQ(ierr);

    // write body forces
    ierr = writeForcesASCII(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // write

// create additional operators (PETSc Mat objects) for the decoupled IBPM
PetscErrorCode DecoupledIBPMSolver::createExtraOperators()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    Mat R;
    Mat MHat;
    Mat BN;
    Vec RDiag;
    Vec MHatDiag;

    // create diagonal matrix R and hold the diagonal in a PETSc Vec object
    ierr = petibm::operators::createR(mesh, R); CHKERRQ(ierr);
    ierr = MatCreateVecs(R, nullptr, &RDiag); CHKERRQ(ierr);
    ierr = MatGetDiagonal(R, RDiag); CHKERRQ(ierr);

    // create diagonal matrix MHat and hold the diagonal in a PETSc Vec object
    ierr = petibm::operators::createMHead(mesh, MHat); CHKERRQ(ierr);
    ierr = MatCreateVecs(MHat, nullptr, &MHatDiag); CHKERRQ(ierr);
    ierr = MatGetDiagonal(MHat, MHatDiag); CHKERRQ(ierr);

    // create a Delta operator and its transpose
    const YAML::Node &node = config["parameters"];
    std::string name = node["delta"].as<std::string>("ROMA_ET_AL_1999");
    petibm::delta::DeltaKernel kernel;
    PetscInt kernelSize;
    ierr = petibm::delta::getKernel(name, kernel, kernelSize); CHKERRQ(ierr);
    ierr = petibm::operators::createDelta(
        mesh, bc, bodies, kernel, kernelSize, E); CHKERRQ(ierr);
    ierr = MatTranspose(E, MAT_INITIAL_MATRIX, &H); CHKERRQ(ierr);

    // create the regularization operator: E
    ierr = MatDiagonalScale(E, nullptr, RDiag); CHKERRQ(ierr);
    ierr = MatDiagonalScale(E, nullptr, MHatDiag); CHKERRQ(ierr);

    // create the spreading operator: H
    // Note: we do nothing here. Because if f is force (unit [m/s^2]), we need
    // to convert it to pressure when using it in Eulerian momentum equations.
    // That is, (H*R^{-1}*f). While H = R*Delta, This means the real scattering
    // operator used in momentum equation is just a Delta operator. So we just
    // let H be the Delta operator. One thing to remember is that this is based
    // on the assumption that the size of Lagrangian element is equal to the
    // size of the Eulerian grid near by.

    // create the operator BN
    PetscInt N;  // order of the truncate Taylor series expansion
    N = config["parameters"]["BN"].as<PetscInt>(1);
    ierr = petibm::operators::createBnHead(
        L, dt, diffCoeffs->implicitCoeff * nu, N, BN); CHKERRQ(ierr);

    // create the operator BNH
    ierr = MatMatMult(
        BN, H, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &BNH); CHKERRQ(ierr);

    // create the operator EBNH
    ierr = MatMatMult(
        E, BNH, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &EBNH); CHKERRQ(ierr);

    // destroy temporary PETSc Vec and Mat objects
    ierr = VecDestroy(&RDiag); CHKERRQ(ierr);
    ierr = VecDestroy(&MHatDiag); CHKERRQ(ierr);
    ierr = MatDestroy(&MHat); CHKERRQ(ierr);
    ierr = MatDestroy(&R); CHKERRQ(ierr);
    ierr = MatDestroy(&BN); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // createExtraOperators

// create additional vectors (PETSc Vec objects) for the decoupled IBPM
PetscErrorCode DecoupledIBPMSolver::createExtraVectors()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = DMCreateGlobalVector(bodies->dmPack, &f); CHKERRQ(ierr);
    ierr = VecDuplicate(f, &df); CHKERRQ(ierr);
    ierr = VecDuplicate(f, &rhsf); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // createExtraVectors

// assemble the right-hand side vector of the forces system
PetscErrorCode DecoupledIBPMSolver::assembleRHSVelocity()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    // assemble the part coming from underlying Navier-Stokes solver
    ierr = NavierStokesSolver::assembleRHSVelocity(); CHKERRQ(ierr);

    ierr = PetscLogStagePush(stageRHSVelocity); CHKERRQ(ierr);

    // add the Lagrangian forces spread to the Eulerian grid
    ierr = MatMultAdd(H, f, rhs1, rhs1); CHKERRQ(ierr);

    ierr = PetscLogStagePop(); CHKERRQ(ierr);  // end of stageRHSVelocity

    PetscFunctionReturn(0);
}  // assembleRHSVelocity

// assemble the RHS vector of the system for the Lagrangian forces
PetscErrorCode DecoupledIBPMSolver::assembleRHSForces()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = PetscLogStagePush(stageRHSForces); CHKERRQ(ierr);

    // rhsf is -E u^{**}
    ierr = MatMult(E, solution->UGlobal, rhsf); CHKERRQ(ierr);
    ierr = VecScale(rhsf, -1.0); CHKERRQ(ierr);

    ierr = PetscLogStagePop(); CHKERRQ(ierr);  // end of stageRHSForces

    PetscFunctionReturn(0);
}  // assembleRHSForces

// solve the system for the increment in the Lagrangian forces
PetscErrorCode DecoupledIBPMSolver::solveForces()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = PetscLogStagePush(stageSolveForces); CHKERRQ(ierr);

    // solve for the increment in the Lagrangian forces
    ierr = fSolver->solve(df, rhsf); CHKERRQ(ierr);

    ierr = PetscLogStagePop(); CHKERRQ(ierr);  // end of stageSolveForces

    PetscFunctionReturn(0);
}  // solveForces

// update the velocity field to satisfy the no-slip condition on the body interfaces
PetscErrorCode DecoupledIBPMSolver::applyNoSlip()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    // u = u + BN H df
    ierr = MatMultAdd(
        BNH, df, solution->UGlobal, solution->UGlobal); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // applyNoSlip

// update the vector forces
PetscErrorCode DecoupledIBPMSolver::updateForces()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = PetscLogStagePush(stageUpdate); CHKERRQ(ierr);

    // f = f + df
    ierr = VecAXPY(f, 1.0, df); CHKERRQ(ierr);

    ierr = PetscLogStagePop(); CHKERRQ(ierr);  // end of stageUpdate

    PetscFunctionReturn(0);
}  // updateForces

// write data required to restart a simulation into a HDF5 file
PetscErrorCode DecoupledIBPMSolver::writeRestartDataHDF5(
    const std::string &filePath)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = NavierStokesSolver::writeRestartDataHDF5(filePath); CHKERRQ(ierr);

    ierr = PetscLogStagePush(stageWrite); CHKERRQ(ierr);

    // create PetscViewer object with append mode
    PetscViewer viewer;
    ierr = PetscViewerCreate(comm, &viewer); CHKERRQ(ierr);
    ierr = PetscViewerSetType(viewer, PETSCVIEWERHDF5); CHKERRQ(ierr);
    ierr = PetscViewerFileSetMode(viewer, FILE_MODE_APPEND); CHKERRQ(ierr);
    ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);

    // go to the root node first (just in case, not necessary)
    ierr = PetscViewerHDF5PushGroup(viewer, "/"); CHKERRQ(ierr);

    // write the Lagrangian forces
    ierr = PetscObjectSetName((PetscObject)f, "force"); CHKERRQ(ierr);
    ierr = VecView(f, viewer); CHKERRQ(ierr);

    // destroy viewer
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

    ierr = PetscLogStagePop(); CHKERRQ(ierr);  // end of stageWrite

    PetscFunctionReturn(0);
}  // writeRestartDataHDF5

// read data required to restart a simulation from a HDF5 file
PetscErrorCode DecoupledIBPMSolver::readRestartDataHDF5(
    const std::string &filePath)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = NavierStokesSolver::readRestartDataHDF5(filePath); CHKERRQ(ierr);

    // create PetscViewer object with read-only mode
    PetscViewer viewer;
    ierr = PetscViewerCreate(comm, &viewer); CHKERRQ(ierr);
    ierr = PetscViewerSetType(viewer, PETSCVIEWERHDF5); CHKERRQ(ierr);
    ierr = PetscViewerFileSetMode(viewer, FILE_MODE_READ); CHKERRQ(ierr);
    ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);

    // go to the root node first (just in case, not necessary)
    ierr = PetscViewerHDF5PushGroup(viewer, "/"); CHKERRQ(ierr);

    // read the Lagrangian forces
    ierr = PetscObjectSetName((PetscObject)f, "force"); CHKERRQ(ierr);
    ierr = VecLoad(f, viewer); CHKERRQ(ierr);

    // destroy viewer
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // readRestartDataHDF5

// write numbers of iterations and residuals of linear solvers to an ASCII file
PetscErrorCode DecoupledIBPMSolver::writeLinSolversInfo()
{
    PetscErrorCode ierr;
    PetscInt nIters;
    PetscReal res;

    PetscFunctionBeginUser;

    ierr = PetscLogStagePush(stageWrite); CHKERRQ(ierr);

    // write the time value
    ierr = PetscViewerASCIIPrintf(solversViewer, "%d\t", ite); CHKERRQ(ierr);

    // write iterations number and residual for the velocity solver
    ierr = vSolver->getIters(nIters); CHKERRQ(ierr);
    ierr = vSolver->getResidual(res); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(
        solversViewer, "%d\t%e\t", nIters, res); CHKERRQ(ierr);

    // write iterations number and residual for the Poisson solver
    ierr = pSolver->getIters(nIters); CHKERRQ(ierr);
    ierr = pSolver->getResidual(res); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(
        solversViewer, "%d\t%e\t", nIters, res); CHKERRQ(ierr);

    // write iterations number and residual for the forces solver
    ierr = fSolver->getIters(nIters); CHKERRQ(ierr);
    ierr = fSolver->getResidual(res); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(
        solversViewer, "%d\t%e\n", nIters, res); CHKERRQ(ierr);

    ierr = PetscLogStagePop(); CHKERRQ(ierr);  // end of stageWrite

    PetscFunctionReturn(0);
}  // writeLinSolversInfo

// integrate the forces and output to ASCII file
PetscErrorCode DecoupledIBPMSolver::writeForcesASCII()
{
    PetscErrorCode ierr;
    petibm::type::RealVec2D fAvg;

    PetscFunctionBeginUser;

    ierr = PetscLogStagePush(stageIntegrateForces); CHKERRQ(ierr);

    // get averaged forces first
    ierr = bodies->calculateAvgForces(f, fAvg); CHKERRQ(ierr);

    ierr = PetscLogStagePop(); CHKERRQ(ierr);

    ierr = PetscLogStagePush(stageWrite); CHKERRQ(ierr);

    // write the time value
    ierr = PetscViewerASCIIPrintf(forcesViewer, "%10.8e\t", t); CHKERRQ(ierr);

    // write forces for each immersed body
    for (int i = 0; i < bodies->nBodies; ++i)
    {
        for (int d = 0; d < mesh->dim; ++d)
        {
            ierr = PetscViewerASCIIPrintf(
                forcesViewer, "%10.8e\t", fAvg[i][d]); CHKERRQ(ierr);
        }
    }
    ierr = PetscViewerASCIIPrintf(forcesViewer, "\n"); CHKERRQ(ierr);

    ierr = PetscLogStagePop(); CHKERRQ(ierr);  // end of stageIntegrateForces

    PetscFunctionReturn(0);
}  // writeForcesASCII
