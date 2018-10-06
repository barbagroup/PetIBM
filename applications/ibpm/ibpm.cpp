/**
 * \file ibpm.cpp
 * \brief Implementation of the class \c IBPMSolver.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 * \see ibpm
 * \ingroup ibpm
 */

#include <petscviewerhdf5.h>

#include <petibm/io.h>

#include "ibpm.h"

IBPMSolver::IBPMSolver(const MPI_Comm &world,
                       const YAML::Node &node)
{
    init(world, node);
}  // IBPMSolver

IBPMSolver::~IBPMSolver()
{
    PetscErrorCode ierr;
    PetscBool finalized;

    PetscFunctionBeginUser;

    ierr = PetscFinalized(&finalized); CHKERRV(ierr);
    if (finalized) return;

    ierr = destroy(); CHKERRV(ierr);
}  // ~IBPMSolver

// manual destroy data
PetscErrorCode IBPMSolver::destroy()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    bodies.reset();
    ierr = ISDestroy(&isDE[0]); CHKERRQ(ierr);
    ierr = ISDestroy(&isDE[1]); CHKERRQ(ierr);
    ierr = VecResetArray(phi); CHKERRQ(ierr);
    ierr = VecDestroy(&phi); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&forcesViewer); CHKERRQ(ierr);
    ierr = NavierStokesSolver::destroy(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // destroy

PetscErrorCode IBPMSolver::init(const MPI_Comm &world, const YAML::Node &node)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = NavierStokesSolver::init(world, node); CHKERRQ(ierr);

    ierr = PetscLogStagePush(stageInitialize); CHKERRQ(ierr);

    // create a pack of immersed bodies
    ierr = petibm::body::createBodyPack(
        comm, mesh->dim, config, bodies); CHKERRQ(ierr);
    ierr = bodies->updateMeshIdx(mesh); CHKERRQ(ierr);

    // create an ASCII PetscViewer to output the body forces
    ierr = createPetscViewerASCII(
        config["directory"].as<std::string>() +
        "/forces-" + std::to_string(ite) + ".txt",
        FILE_MODE_WRITE, forcesViewer); CHKERRQ(ierr);

    // register additional logging stage
    ierr = PetscLogStageRegister(
        "integrateForces", &stageIntegrateForces); CHKERRQ(ierr);

    ierr = PetscLogStagePop(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // init

// write solution fields, linear solvers info, and body forces to files
PetscErrorCode IBPMSolver::write()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = NavierStokesSolver::write(); CHKERRQ(ierr);

    // write body forces
    ierr = writeForcesASCII(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // write

// create the linear operators of the solver (PETSc Mat objects)
PetscErrorCode IBPMSolver::createOperators()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    Mat R, MHat, BN, GH[2], DE[2];  // temporary operators
    Vec RDiag, MHatDiag;            // temporary vectors
    IS is[2];                       // temporary index sets

    // create the divergence operator: D
    ierr = petibm::operators::createDivergence(
        mesh, bc, DE[0], DCorrection, PETSC_FALSE); CHKERRQ(ierr);
    
    // create the gradient operator: G
    ierr = petibm::operators::createGradient(
        mesh, GH[0], PETSC_FALSE); CHKERRQ(ierr);
    
    // create the Laplacian operator: L
    ierr = petibm::operators::createLaplacian(
        mesh, bc, L, LCorrection); CHKERRQ(ierr);
    
    // create the operator for the convective terms: N
    ierr = petibm::operators::createConvection(
        mesh, bc, N); CHKERRQ(ierr);

    // create the implicit operator for the velocity system: A
    ierr = MatDuplicate(L, MAT_COPY_VALUES, &A); CHKERRQ(ierr);
    ierr = MatScale(A, -diffCoeffs->implicitCoeff * nu); CHKERRQ(ierr);
    ierr = MatShift(A, 1.0 / dt); CHKERRQ(ierr);

    // create diagonal matrix R and hold the diagonal in a PETSc Vec object
    ierr = petibm::operators::createR(mesh, R); CHKERRQ(ierr);
    ierr = MatCreateVecs(R, nullptr, &RDiag); CHKERRQ(ierr);
    ierr = MatGetDiagonal(R, RDiag); CHKERRQ(ierr);

    // create diagonal matrix MHat and hold the diagonal in a PETSc Vec object
    ierr = petibm::operators::createMHead(mesh, MHat); CHKERRQ(ierr);
    ierr = MatCreateVecs(MHat, nullptr, &MHatDiag); CHKERRQ(ierr);
    ierr = MatGetDiagonal(MHat, MHatDiag); CHKERRQ(ierr);

    // create a Delta operator and its transpose (equal to H)
    ierr = petibm::operators::createDelta(
        mesh, bc, bodies, DE[1]); CHKERRQ(ierr);
    ierr = MatTranspose(DE[1], MAT_INITIAL_MATRIX, &GH[1]); CHKERRQ(ierr);

    // create the regularization operator: E
    ierr = MatDiagonalScale(DE[1], nullptr, RDiag); CHKERRQ(ierr);
    ierr = MatDiagonalScale(DE[1], nullptr, MHatDiag); CHKERRQ(ierr);
    ierr = VecDestroy(&RDiag); CHKERRQ(ierr);
    ierr = VecDestroy(&MHatDiag); CHKERRQ(ierr);
    ierr = MatDestroy(&R); CHKERRQ(ierr);
    ierr = MatDestroy(&MHat); CHKERRQ(ierr);

    // create the opposite of the spreading operator: -H
    ierr = MatScale(GH[1], -1.0); CHKERRQ(ierr);

    // get combined operator G; R is used temporarily
    ierr = MatCreateNest(
        comm, 1, nullptr, 2, nullptr, GH, &R); CHKERRQ(ierr);
    ierr = MatConvert(R, MATAIJ, MAT_INITIAL_MATRIX, &G); CHKERRQ(ierr);
    ierr = MatDestroy(&R); CHKERRQ(ierr);
    ierr = MatDestroy(&GH[0]); CHKERRQ(ierr);
    ierr = MatDestroy(&GH[1]); CHKERRQ(ierr);

    // get combined operator D; R is used temporarily; also get ISs
    ierr = MatCreateNest(
        comm, 2, nullptr, 1, nullptr, DE, &R); CHKERRQ(ierr);
    ierr = MatConvert(R, MATAIJ, MAT_INITIAL_MATRIX, &D); CHKERRQ(ierr);
    ierr = MatNestGetISs(R, is, nullptr); CHKERRQ(ierr);
    ierr = ISDuplicate(is[0], &isDE[0]); CHKERRQ(ierr);
    ierr = ISDuplicate(is[1], &isDE[1]); CHKERRQ(ierr);
    ierr = ISCopy(is[0], isDE[0]); CHKERRQ(ierr);
    ierr = ISCopy(is[1], isDE[1]); CHKERRQ(ierr);
    ierr = MatDestroy(&R); CHKERRQ(ierr);
    ierr = MatDestroy(&DE[0]); CHKERRQ(ierr);
    ierr = MatDestroy(&DE[1]); CHKERRQ(ierr);

    // create the projection operator: BNG
    ierr = petibm::operators::createBnHead(
        L, dt, diffCoeffs->implicitCoeff * nu, 1, BN); CHKERRQ(ierr);
    ierr = MatMatMult(
        BN, G, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &BNG); CHKERRQ(ierr);

    // create the modified Poisson operator
    ierr = MatMatMult(
        D, BNG, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &DBNG); CHKERRQ(ierr);

    // set the nullspace of the modified Poisson system
    ierr = setNullSpace(); CHKERRQ(ierr);

    // destroy temporary operator
    ierr = MatDestroy(&BN); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // createOperators

// create the vectors of the solver (PETSc Vec objects)
PetscErrorCode IBPMSolver::createVectors()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    Vec temp;
    const PetscReal *data;

    // create the vector of the couple (pressure field, Lagrangian forces)
    ierr = MatCreateVecs(G, &phi, nullptr); CHKERRQ(ierr);

    // swap pGlobal and phi to reuse functions from the Navier-Stokes solver
    temp = solution->pGlobal;
    solution->pGlobal = phi;
    phi = temp;
    temp = PETSC_NULL;

    // destroy phi's underlying raw array but keep all other information
    ierr = VecReplaceArray(phi, nullptr); CHKERRQ(ierr);

    // reset the underlying data of phi to the pressure portion in pGlobal
    ierr = VecGetSubVector(solution->pGlobal, isDE[0], &temp); CHKERRQ(ierr);
    ierr = VecGetArrayRead(temp, &data); CHKERRQ(ierr);
    ierr = VecPlaceArray(phi, data); CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(temp, &data); CHKERRQ(ierr);
    ierr = VecRestoreSubVector(
        solution->pGlobal, isDE[0], &temp); CHKERRQ(ierr);

    // create remaining vectors
    ierr = NavierStokesSolver::createVectors(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // createVectors

// set the nullspace of the modified Poisson system
PetscErrorCode IBPMSolver::setNullSpace()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    std::string type;
    ierr = pSolver->getType(type); CHKERRQ(ierr);

    if (type == "PETSc KSP")
    {
        Vec n, phiPortion;
        MatNullSpace nsp;
        ierr = MatCreateVecs(DBNG, &n, nullptr); CHKERRQ(ierr);
        ierr = VecSet(n, 0.0); CHKERRQ(ierr);
        ierr = VecGetSubVector(n, isDE[0], &phiPortion); CHKERRQ(ierr);
        ierr = VecSet(phiPortion, 1.0 / std::sqrt(mesh->pN)); CHKERRQ(ierr);
        ierr = VecRestoreSubVector(n, isDE[0], &phiPortion); CHKERRQ(ierr);
        ierr = MatNullSpaceCreate(
            comm, PETSC_FALSE, 1, &n, &nsp); CHKERRQ(ierr);
        ierr = MatSetNullSpace(DBNG, nsp); CHKERRQ(ierr);
        ierr = MatSetNearNullSpace(DBNG, nsp); CHKERRQ(ierr);
        ierr = VecDestroy(&n); CHKERRQ(ierr);
        ierr = MatNullSpaceDestroy(&nsp); CHKERRQ(ierr);
        isRefP = PETSC_FALSE;
    }
    else if (type == "NVIDIA AmgX")
    {
        PetscInt row[1] = {0};
        ierr = MatZeroRowsColumns(
            DBNG, 1, row, 1.0, nullptr, nullptr); CHKERRQ(ierr);
        isRefP = PETSC_TRUE;
    }
    else
    {
        SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                 "Could not recognize the type of linear solver: %s\n",
                 type.c_str());
    }

    PetscFunctionReturn(0);
}  // setNullSpace

// assemble the right-hand side vector of the modified Poisson system
PetscErrorCode IBPMSolver::assembleRHSPoisson()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = PetscLogStagePush(stageRHSPoisson); CHKERRQ(ierr);

    // compute the divergence of the intermediate velocity field
    ierr = MatMult(D, solution->UGlobal, rhs2); CHKERRQ(ierr);

    // note: bc2 is a subset of rhs2
    Vec bc2;
    ierr = VecGetSubVector(rhs2, isDE[0], &bc2); CHKERRQ(ierr);
    ierr = MatMultAdd(DCorrection, solution->UGlobal, bc2, bc2); CHKERRQ(ierr);
    ierr = VecRestoreSubVector(rhs2, isDE[0], &bc2); CHKERRQ(ierr);

    if (isRefP)  // if the pressure is pinned at one point
    {
        ierr = VecSetValue(rhs2, 0, 0.0, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecAssemblyBegin(rhs2); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(rhs2); CHKERRQ(ierr);
    }

    ierr = PetscLogStagePop(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // assembleRHSPoisson

// write the solution fields into a HDF5 file
PetscErrorCode IBPMSolver::writeSolutionHDF5(const std::string &filePath)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    Vec temp = PETSC_NULL;

    // let solution->pGlobal point to phi, so that we can use solution->write
    temp = solution->pGlobal;
    solution->pGlobal = phi;

    ierr = NavierStokesSolver::writeSolutionHDF5(filePath); CHKERRQ(ierr);

    // restore pointers
    solution->pGlobal = temp;
    temp = PETSC_NULL;

    PetscFunctionReturn(0);
}  // writeSolutionHDF5

// write all data required to restart a simulation into a HDF5 file
PetscErrorCode IBPMSolver::writeRestartDataHDF5(const std::string &filePath)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = NavierStokesSolver::writeRestartDataHDF5(filePath); CHKERRQ(ierr);

    // write forces
    Vec f;
    ierr = VecGetSubVector(solution->pGlobal, isDE[1], &f); CHKERRQ(ierr);
    ierr = petibm::io::writeHDF5Vecs(
        comm, filePath, "/", {"force"}, {f}, FILE_MODE_APPEND); CHKERRQ(ierr);
    ierr = VecRestoreSubVector(solution->pGlobal, isDE[1], &f); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // writeRestartDataHDF5

// read all data required to restart a simulation into a HDF5 file
PetscErrorCode IBPMSolver::readRestartDataHDF5(const std::string &filePath)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    Vec temp = PETSC_NULL;

    // let solution->pGlobal point to phi, so that we can use solution->write
    temp = solution->pGlobal;
    solution->pGlobal = phi;

    ierr = NavierStokesSolver::readRestartDataHDF5(filePath); CHKERRQ(ierr);

    // restore pointers
    solution->pGlobal = temp;
    temp = PETSC_NULL;

    // write forces
    std::vector<Vec> f(1);
    ierr = VecGetSubVector(solution->pGlobal, isDE[1], &f[0]); CHKERRQ(ierr);
    ierr = petibm::io::readHDF5Vecs(
        comm, filePath, "/", {"force"}, f); CHKERRQ(ierr);
    ierr = VecRestoreSubVector(
        solution->pGlobal, isDE[1], &f[0]); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // readRestartDataHDF5

// integrate the forces and output to ASCII file
PetscErrorCode IBPMSolver::writeForcesASCII()
{
    PetscErrorCode ierr;
    petibm::type::RealVec2D fAvg;

    PetscFunctionBeginUser;

    ierr = PetscLogStagePush(stageIntegrateForces); CHKERRQ(ierr);

    // get sub section f and calculate averaged forces
    Vec f;
    ierr = VecGetSubVector(solution->pGlobal, isDE[1], &f); CHKERRQ(ierr);
    ierr = bodies->calculateAvgForces(f, fAvg); CHKERRQ(ierr);
    ierr = VecRestoreSubVector(solution->pGlobal, isDE[1], &f); CHKERRQ(ierr);

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

    ierr = PetscLogStagePop(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // writeForcesASCII
