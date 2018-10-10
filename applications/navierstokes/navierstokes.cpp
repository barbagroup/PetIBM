/**
 * \file navierstokes.cpp
 * \brief Implementation of the class \c NavierStokesSolver.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 * \see nssolver
 * \ingroup nssolver
 */

#include <iomanip>

#include <petscviewerhdf5.h>

#include <petibm/io.h>

#include "navierstokes.h"

NavierStokesSolver::NavierStokesSolver(const MPI_Comm &world,
                                       const YAML::Node &node)
{
    init(world, node);
}  // NavierStokesSolver

NavierStokesSolver::~NavierStokesSolver()
{
    PetscErrorCode ierr;
    PetscBool finalized;

    PetscFunctionBeginUser;

    ierr = PetscFinalized(&finalized); CHKERRV(ierr);
    if (finalized) return;

    ierr = destroy(); CHKERRV(ierr);
}  // ~NavierStokesSolver

// destroy objects of the NavierStokesSolver class
PetscErrorCode NavierStokesSolver::destroy()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    comm = MPI_COMM_NULL;
    commSize = commRank = 0;

    // destroy vectors of the solver (PETSc Vec objects)
    ierr = VecDestroy(&dP); CHKERRQ(ierr);
    ierr = VecDestroy(&bc1); CHKERRQ(ierr);
    ierr = VecDestroy(&rhs1); CHKERRQ(ierr);
    ierr = VecDestroy(&rhs2); CHKERRQ(ierr);
    for (unsigned int i = 0; i < conv.size(); ++i)
    {
        ierr = VecDestroy(&conv[i]); CHKERRQ(ierr);
    }
    for (unsigned int i = 0; i < diff.size(); ++i)
    {
        ierr = VecDestroy(&diff[i]); CHKERRQ(ierr);
    }

    // destroy operators of the solver (PETSc Mat objects)
    ierr = MatDestroy(&A); CHKERRQ(ierr);
    ierr = MatDestroy(&DBNG); CHKERRQ(ierr);
    ierr = MatDestroy(&BNG); CHKERRQ(ierr);
    ierr = MatDestroy(&N); CHKERRQ(ierr);
    ierr = MatDestroy(&G); CHKERRQ(ierr);
    ierr = MatDestroy(&D); CHKERRQ(ierr);
    ierr = MatDestroy(&DCorrection); CHKERRQ(ierr);
    ierr = MatDestroy(&L); CHKERRQ(ierr);
    ierr = MatDestroy(&LCorrection); CHKERRQ(ierr);

    // destroy the probes
    for (auto probe : probes)
    {
        ierr = probe->destroy(); CHKERRQ(ierr);
    }

    // decrease reference count or destroy
    vSolver.reset();
    pSolver.reset();
    config.reset();
    convCoeffs.reset();
    diffCoeffs.reset();
    solution.reset();
    bc.reset();
    mesh.reset();

    ierr = PetscViewerDestroy(&solversViewer); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // destroy

PetscErrorCode NavierStokesSolver::init(const MPI_Comm &world, const YAML::Node &node)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = PetscLogStageRegister(
        "initialize", &stageInitialize); CHKERRQ(ierr);
    ierr = PetscLogStagePush(stageInitialize); CHKERRQ(ierr);

    // record the MPI communicator, size, and process rank
    comm = world;
    ierr = MPI_Comm_size(comm, &commSize); CHKERRQ(ierr);
    ierr = MPI_Comm_rank(comm, &commRank); CHKERRQ(ierr);

    // record the YAML configuration settings
    config = node;

    // get the time-step size
    dt = config["parameters"]["dt"].as<PetscReal>();
    // get the starting time step and time value
    nstart = config["parameters"]["startStep"].as<PetscInt>(0);
    ite = nstart;
    t = config["parameters"]["t"].as<PetscReal>(0.0);
    // get the number of time steps to compute
    nt = config["parameters"]["nt"].as<PetscInt>();
    // get the saving frequencies
    nsave = config["parameters"]["nsave"].as<PetscInt>();
    nrestart = config["parameters"]["nrestart"].as<PetscInt>();
    // get the viscous diffusion coefficient
    nu = config["flow"]["nu"].as<PetscReal>();

    // create the Cartesian mesh
    ierr = petibm::mesh::createMesh(comm, config, mesh); CHKERRQ(ierr);
    // write the grid points into a HDF5 file
    std::string filePath = config["directory"].as<std::string>() + "/grid.h5";
    ierr = mesh->write(filePath); CHKERRQ(ierr);

    // create the data object for the boundary conditions
    ierr = petibm::boundary::createBoundary(mesh, config, bc); CHKERRQ(ierr);

    // create the solution object
    ierr = petibm::solution::createSolution(mesh, solution); CHKERRQ(ierr);

    // set the initial conditions
    ierr = solution->setInitialConditions(config); CHKERRQ(ierr);

    // initialize ghost-point values and equations;
    // must be done before creating the operators
    ierr = bc->setGhostICs(solution); CHKERRQ(ierr);

    // create the time-scheme objects
    ierr = petibm::timeintegration::createTimeIntegration(
        "convection", config, convCoeffs); CHKERRQ(ierr);
    ierr = petibm::timeintegration::createTimeIntegration(
        "diffusion", config, diffCoeffs); CHKERRQ(ierr);

    // create the linear solver objects
    ierr = petibm::linsolver::createLinSolver(
        "velocity", config, vSolver); CHKERRQ(ierr);
    ierr = petibm::linsolver::createLinSolver(
        "poisson", config, pSolver); CHKERRQ(ierr);

    // create operators (PETSc Mat objects)
    ierr = createOperators(); CHKERRQ(ierr);

    // create PETSc Vec objects
    ierr = createVectors(); CHKERRQ(ierr);

    // set coefficient matrix of the linear solvers
    ierr = vSolver->setMatrix(A); CHKERRQ(ierr);
    ierr = pSolver->setMatrix(DBNG); CHKERRQ(ierr);

    // create probes to monitor the solution in some regions of the domain
    probes.resize(config["probes"].size());
    for (unsigned int i = 0; i < probes.size(); ++i)
    {
        ierr = petibm::misc::createProbe(mesh->comm, config["probes"][i],
                                         mesh, probes[i]); CHKERRQ(ierr);
    }

    // create an ASCII PetscViewer to output linear solvers info
    ierr = createPetscViewerASCII(
        config["directory"].as<std::string>() +
        "/iterations-" + std::to_string(ite) + ".txt",
        FILE_MODE_WRITE, solversViewer); CHKERRQ(ierr);

    // register logging stages
    ierr = PetscLogStageRegister(
        "rhsVelocity", &stageRHSVelocity); CHKERRQ(ierr);
    ierr = PetscLogStageRegister(
        "solveVelocity", &stageSolveVelocity); CHKERRQ(ierr);
    ierr = PetscLogStageRegister(
        "rhsPoisson", &stageRHSPoisson); CHKERRQ(ierr);
    ierr = PetscLogStageRegister(
        "solvePoisson", &stageSolvePoisson); CHKERRQ(ierr);
    ierr = PetscLogStageRegister(
        "projectionStep", &stageProjectionStep); CHKERRQ(ierr);
    ierr = PetscLogStageRegister(
        "write", &stageWrite); CHKERRQ(ierr);
    ierr = PetscLogStageRegister(
        "monitor", &stageMonitor); CHKERRQ(ierr);

    ierr = PetscLogStagePop(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // init

// read or write initial data
PetscErrorCode NavierStokesSolver::ioInitialData()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    if (ite == 0)  // write the initial solution fields to a HDF5 file
    {
        std::stringstream ss;
        std::string filePath;
        ss << std::setfill('0') << std::setw(7) << ite;
        filePath = config["output"].as<std::string>() + "/" + ss.str() + ".h5";
        ierr = PetscPrintf(comm, "[time step %d] Writing solution data... ",
                            ite); CHKERRQ(ierr);
        ierr = writeSolutionHDF5(filePath); CHKERRQ(ierr);
        ierr = PetscPrintf(comm, "done\n"); CHKERRQ(ierr);
    }
    else  // read restart data from HDF5 file
    {
        std::stringstream ss;
        std::string filePath;
        ss << std::setfill('0') << std::setw(7) << ite;
        filePath = config["output"].as<std::string>() + "/" + ss.str() + ".h5";
        ierr = PetscPrintf(comm, "[time step %d] Reading restart data... ",
                            ite); CHKERRQ(ierr);
        ierr = readRestartDataHDF5(filePath); CHKERRQ(ierr);
        ierr = PetscPrintf(comm, "done\n"); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}  // ioInitialData

// advance the flow solver by one time-step
PetscErrorCode NavierStokesSolver::advance()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    t += dt;
    ite++;

    // prepare velocity system and solve it
    ierr = assembleRHSVelocity(); CHKERRQ(ierr);
    ierr = solveVelocity(); CHKERRQ(ierr);

    // prepare Poisson system and solve it
    ierr = assembleRHSPoisson(); CHKERRQ(ierr);
    ierr = solvePoisson(); CHKERRQ(ierr);

    // project the velocity field onto the divergence-free space
    ierr = projectionStep(); CHKERRQ(ierr);

    // update ghost-point values
    ierr = bc->updateGhostValues(solution); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // advance

// write solutions fields and linear solvers info to files
PetscErrorCode NavierStokesSolver::write()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    // write linear solvers info
    ierr = writeLinSolversInfo(); CHKERRQ(ierr);

    if (ite % nsave == 0)  // write solution fields
    {
        std::stringstream ss;
        std::string filePath;
        ss << std::setfill('0') << std::setw(7) << ite;
        filePath = config["output"].as<std::string>() + "/" + ss.str() + ".h5";
        ierr = PetscPrintf(comm, "[time step %d] Writing solution data... ",
                           ite); CHKERRQ(ierr);
        ierr = writeSolutionHDF5(filePath); CHKERRQ(ierr);
        ierr = PetscPrintf(comm, "done\n"); CHKERRQ(ierr);
        // output the PETSc log to an ASCII file
        filePath = config["logs"].as<std::string>() + "/" + ss.str() + ".log";
        ierr = petibm::io::writePetscLog(comm, filePath); CHKERRQ(ierr);
    }
    if (ite % nrestart == 0)  // write restart data
    {
        std::string filePath;
        std::stringstream ss;
        ss << std::setfill('0') << std::setw(7) << ite;
        filePath = config["output"].as<std::string>() + "/" + ss.str() + ".h5";
        ierr = PetscPrintf(comm, "[time step %d] Writing restart data... ",
                           ite); CHKERRQ(ierr);
        ierr = writeRestartDataHDF5(filePath); CHKERRQ(ierr);
        ierr = PetscPrintf(comm, "done\n"); CHKERRQ(ierr);
    }

    // monitor probes and write to files
    ierr = monitorProbes(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // write

// evaluate if the simulation is finished
bool NavierStokesSolver::finished()
{
    return ite >= nstart + nt;
}  // finished

// create the linear operators of the solver (PETSc Mat objects)
PetscErrorCode NavierStokesSolver::createOperators()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    Mat BN;  // a temporary operator

    // create the divergence operator: D
    ierr = petibm::operators::createDivergence(
        mesh, bc, D, DCorrection, PETSC_FALSE); CHKERRQ(ierr);
    
    // create the gradient operator: G
    ierr = petibm::operators::createGradient(
        mesh, G, PETSC_FALSE); CHKERRQ(ierr);
    
    // create the Laplacian operator: L
    ierr = petibm::operators::createLaplacian(
        mesh, bc, L, LCorrection); CHKERRQ(ierr);
    
    // create the operator for the convective terms: N
    ierr = petibm::operators::createConvection(
        mesh, bc, N); CHKERRQ(ierr);

    // create the implicit operator of the velocity system: A
    ierr = MatDuplicate(L, MAT_COPY_VALUES, &A); CHKERRQ(ierr);
    ierr = MatScale(A, -diffCoeffs->implicitCoeff * nu); CHKERRQ(ierr);
    ierr = MatShift(A, 1.0 / dt); CHKERRQ(ierr);

    // create the projection operator: BNG
    ierr = petibm::operators::createBnHead(
        L, dt, diffCoeffs->implicitCoeff * nu, 1, BN); CHKERRQ(ierr);
    ierr = MatMatMult(
        BN, G, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &BNG); CHKERRQ(ierr);

    // create the Poisson operator: DBNG
    ierr = MatMatMult(
        D, BNG, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &DBNG); CHKERRQ(ierr);

    // set the nullspace of the Poisson system
    ierr = setNullSpace(); CHKERRQ(ierr);

    // destroy the temporary operator
    ierr = MatDestroy(&BN); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // createOperators

// create the vectors of the solver (PETSc Vec objects)
PetscErrorCode NavierStokesSolver::createVectors()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = VecDuplicate(solution->pGlobal, &dP); CHKERRQ(ierr);
    ierr = VecDuplicate(solution->UGlobal, &bc1); CHKERRQ(ierr);
    ierr = VecDuplicate(solution->UGlobal, &rhs1); CHKERRQ(ierr);
    ierr = VecDuplicate(solution->pGlobal, &rhs2); CHKERRQ(ierr);

    conv.resize(convCoeffs->nExplicit);
    for (unsigned int i = 0; i < conv.size(); ++i)
    {
        ierr = VecDuplicate(solution->UGlobal, &conv[i]); CHKERRQ(ierr);
    }

    diff.resize(diffCoeffs->nExplicit);
    for (unsigned int i = 0; i < diff.size(); ++i)
    {
        ierr = VecDuplicate(solution->UGlobal, &diff[i]); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}  // createVectors

// set the nullspace of the Poisson system
PetscErrorCode NavierStokesSolver::setNullSpace()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    std::string type;
    ierr = pSolver->getType(type); CHKERRQ(ierr);

    if (type == "PETSc KSP")
    {
        MatNullSpace nsp;
        ierr = MatNullSpaceCreate(
            comm, PETSC_TRUE, 0, nullptr, &nsp); CHKERRQ(ierr);
        ierr = MatSetNullSpace(DBNG, nsp); CHKERRQ(ierr);
        ierr = MatSetNearNullSpace(DBNG, nsp); CHKERRQ(ierr);
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

// assemble the right-hand side vector of the velocity system
PetscErrorCode NavierStokesSolver::assembleRHSVelocity()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = PetscLogStagePush(stageRHSVelocity); CHKERRQ(ierr);

    // initialize RHS vector with pressure gradient at time-step n
    // $rhs_1 = - \frac{\partial p^n}{\partial x}$
    ierr = MatMult(G, solution->pGlobal, rhs1); CHKERRQ(ierr);
    ierr = VecScale(rhs1, -1.0); CHKERRQ(ierr);

    // add explicit part of time derivative to the RHS vector
    // $rhs_1 += \frac{u^n}{\Delta t}$
    ierr = VecAXPY(rhs1, 1.0 / dt, solution->UGlobal); CHKERRQ(ierr);

    // add all explicit convective terms to the RHS vector
    // $rhs_1 += \sum_{k=0}^s conv_{n - k}$
    {
        // 1. discard the term at the oldest time-step
        // and decrease the time-step by 1
        for (int i = conv.size() - 1; i > 0; i--)
        {
            ierr = VecSwap(conv[i], conv[i - 1]); CHKERRQ(ierr);
        }

        if (conv.size() > 0)
        {
            // 2. compute and store the convective term from time-step n
            ierr = MatMult(N, solution->UGlobal, conv[0]); CHKERRQ(ierr);

            // 3. scale the newest convective term by -1
            // (to be added to the RHS vector)
            ierr = VecScale(conv[0], -1.0); CHKERRQ(ierr);
        }

        // 4. add all explicit convective terms to the RHS vector
        for (unsigned int i = 0; i < conv.size(); ++i)
        {
            ierr = VecAXPY(rhs1, convCoeffs->explicitCoeffs[i], conv[i]);
            CHKERRQ(ierr);
        }
    }

    // add all explicit diffusion terms to the RHS vector
    // $rhs_1 += \sum_{k=0}^s diff_{n - k}$
    {
        // 1. discard the term at the oldest time-step
        // and decrease the time-step by 1
        for (int i = diff.size() - 1; i > 0; i--)
        {
            ierr = VecSwap(diff[i], diff[i - 1]); CHKERRQ(ierr);
        }

        // 2. compute and store the diffusion term from time-step n
        if (diff.size() > 0)
        {
            ierr = MatMult(L, solution->UGlobal, diff[0]); CHKERRQ(ierr);
            ierr = MatMultAdd(
                LCorrection, solution->UGlobal, diff[0], diff[0]); CHKERRQ(ierr);
            ierr = VecScale(diff[0], nu); CHKERRQ(ierr);
        }

        // 3. add all explicit diffusion terms to the RHS vector
        for (unsigned int i = 0; i < diff.size(); ++i)
        {
            ierr = VecAXPY(
                rhs1, diffCoeffs->explicitCoeffs[i], diff[i]); CHKERRQ(ierr);
        }
    }

    // add implicit BC correction terms arising from the diffusion term
    // to the RHS vector
    {
        // 1. update the ghost-point equations
        ierr = bc->updateEqs(solution, dt); CHKERRQ(ierr);

        // 2. compute the implicit BC correction terms
        ierr = MatMult(LCorrection, solution->UGlobal, bc1); CHKERRQ(ierr);
        ierr = VecScale(bc1, nu); CHKERRQ(ierr);

        // 3. add the correction terms to the RHS vector
        ierr = VecAXPY(rhs1, diffCoeffs->implicitCoeff, bc1); CHKERRQ(ierr);
    }

    ierr = PetscLogStagePop(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // assembleRHSVelocity

// solve the linear system for the velocity
PetscErrorCode NavierStokesSolver::solveVelocity()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = PetscLogStagePush(stageSolveVelocity); CHKERRQ(ierr);

    ierr = vSolver->solve(solution->UGlobal, rhs1); CHKERRQ(ierr);

    ierr = PetscLogStagePop(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // solveVelocity

// assemble the right-hand side vector of the Poisson system
PetscErrorCode NavierStokesSolver::assembleRHSPoisson()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = PetscLogStagePush(stageRHSPoisson); CHKERRQ(ierr);

    // compute the divergence of the intermediate velocity field
    ierr = MatMult(D, solution->UGlobal, rhs2); CHKERRQ(ierr);
    ierr = MatMultAdd(DCorrection, solution->UGlobal, rhs2, rhs2);
    CHKERRQ(ierr);

    if (isRefP)  // if the pressure is pinned at one point
    {
        ierr = VecSetValue(rhs2, 0, 0.0, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecAssemblyBegin(rhs2); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(rhs2); CHKERRQ(ierr);
    }

    ierr = PetscLogStagePop(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // assembleRHSPoisson

// solve the Poisson linear system for the pressure correction
PetscErrorCode NavierStokesSolver::solvePoisson()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = PetscLogStagePush(stageSolvePoisson); CHKERRQ(ierr);

    // solve for the pressure correction
    ierr = pSolver->solve(dP, rhs2); CHKERRQ(ierr);

    ierr = PetscLogStagePop(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // solvePoisson

// project the velocity field onto the divergence-free space
// and update the pressure field
PetscErrorCode NavierStokesSolver::projectionStep()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = PetscLogStagePush(stageProjectionStep); CHKERRQ(ierr);

    // project the velocity field onto divergence-free space
    ierr = MatMult(BNG, dP, rhs1); CHKERRQ(ierr);
    ierr = VecAXPY(solution->UGlobal, -1.0, rhs1); CHKERRQ(ierr);

    // update the pressure field
    ierr = VecAXPY(solution->pGlobal, 1.0, dP); CHKERRQ(ierr);

    ierr = PetscLogStagePop(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // projectionStep

// write the solution fields into a HDF5 file
PetscErrorCode NavierStokesSolver::writeSolutionHDF5(const std::string &filePath)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = PetscLogStagePush(stageWrite); CHKERRQ(ierr);

    // write the solution fields to a file
    ierr = solution->write(filePath); CHKERRQ(ierr);
    // write the time value as an attribute of the pressure field dataset
    ierr = writeTimeHDF5(t, filePath); CHKERRQ(ierr);

    ierr = PetscLogStagePop(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // writeSolutionHDF5

// write all data required to restart a simulation into a HDF5 file
PetscErrorCode NavierStokesSolver::writeRestartDataHDF5(
    const std::string &filePath)
{
    PetscErrorCode ierr;
    PetscViewer viewer;
    PetscBool fileExist = PETSC_FALSE;

    PetscFunctionBeginUser;

    ierr = PetscLogStagePush(stageWrite); CHKERRQ(ierr);

    // check if file exist
    ierr = PetscTestFile(filePath.c_str(), 'w', &fileExist); CHKERRQ(ierr);

    if (!fileExist)  // if not, create one and write field solution into it
    {
        ierr = writeSolutionHDF5(filePath); CHKERRQ(ierr);
    }

    // create PetscViewer object with append mode
    ierr = PetscViewerCreate(comm, &viewer); CHKERRQ(ierr);
    ierr = PetscViewerSetType(viewer, PETSCVIEWERHDF5); CHKERRQ(ierr);
    ierr = PetscViewerFileSetMode(viewer, FILE_MODE_APPEND); CHKERRQ(ierr);
    ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);

    // go to the root node first (just in case, not necessary)
    ierr = PetscViewerHDF5PushGroup(viewer, "/"); CHKERRQ(ierr);

    // write explicit convective terms
    ierr = PetscViewerHDF5PushGroup(viewer, "/convection"); CHKERRQ(ierr);
    for (unsigned int i = 0; i < conv.size(); ++i)
    {
        ierr = PetscObjectSetName(
            (PetscObject)conv[i], std::to_string(i).c_str()); CHKERRQ(ierr);
        ierr = VecView(conv[i], viewer); CHKERRQ(ierr);
    }
    // write explicit diffusion terms
    ierr = PetscViewerHDF5PushGroup(viewer, "/diffusion"); CHKERRQ(ierr);
    for (unsigned int i = 0; i < diff.size(); ++i)
    {
        ierr = PetscObjectSetName(
            (PetscObject)diff[i], std::to_string(i).c_str()); CHKERRQ(ierr);
        ierr = VecView(diff[i], viewer); CHKERRQ(ierr);
    }

    // destroy viewer
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

    ierr = PetscLogStagePop(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // writeRestartDataHDF5

// read all data required to restart a simulation into a HDF5 file
PetscErrorCode NavierStokesSolver::readRestartDataHDF5(
    const std::string &filePath)
{
    PetscErrorCode ierr;
    PetscViewer viewer;
    PetscBool fileExist = PETSC_FALSE;

    PetscFunctionBeginUser;

    // check if file exist
    ierr = PetscTestFile(filePath.c_str(), 'r', &fileExist); CHKERRQ(ierr);

    if (!fileExist)  // if the  not, return error
        SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_FILE_READ,
                 "Could not find file \"%s\" for restarting.\n",
                 filePath.c_str());

    // read primary fields
    ierr = solution->read(filePath); CHKERRQ(ierr);
    ierr = readTimeHDF5(filePath, t); CHKERRQ(ierr);

    // create a PetscViewer object with append mode
    ierr = PetscViewerCreate(comm, &viewer); CHKERRQ(ierr);
    ierr = PetscViewerSetType(viewer, PETSCVIEWERHDF5); CHKERRQ(ierr);
    ierr = PetscViewerFileSetMode(viewer, FILE_MODE_READ); CHKERRQ(ierr);
    ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);

    // go to the root node first (just in case, not necessary)
    ierr = PetscViewerHDF5PushGroup(viewer, "/"); CHKERRQ(ierr);

    // read explicit convective terms
    ierr = PetscViewerHDF5PushGroup(viewer, "/convection"); CHKERRQ(ierr);
    for (unsigned int i = 0; i < conv.size(); ++i)
    {
        ierr = PetscObjectSetName(
            (PetscObject)conv[i], std::to_string(i).c_str()); CHKERRQ(ierr);
        ierr = VecLoad(conv[i], viewer); CHKERRQ(ierr);
    }
    // read explicit diffusion terms
    ierr = PetscViewerHDF5PushGroup(viewer, "/diffusion"); CHKERRQ(ierr);
    for (unsigned int i = 0; i < diff.size(); ++i)
    {
        ierr = PetscObjectSetName(
            (PetscObject)diff[i], std::to_string(i).c_str()); CHKERRQ(ierr);
        ierr = VecLoad(diff[i], viewer); CHKERRQ(ierr);
    }

    // destroy viewer
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

    // update ghost-point values and equations based on the current solutio
    // TODO: for convective BCs, it's not totally correct
    ierr = bc->setGhostICs(solution); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // readRestartDataHDF5

// initialize an ASCII PetscViewer object
PetscErrorCode NavierStokesSolver::createPetscViewerASCII(
    const std::string &filePath, const PetscFileMode &mode,
    PetscViewer &viewer)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = PetscViewerCreate(comm, &viewer); CHKERRQ(ierr);
    ierr = PetscViewerSetType(viewer, PETSCVIEWERASCII); CHKERRQ(ierr);
    ierr = PetscViewerFileSetMode(viewer, mode); CHKERRQ(ierr);
    ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // createPetscViewerASCII

// write numbers of iterations and residuals of linear solvers to an ASCII file
PetscErrorCode NavierStokesSolver::writeLinSolversInfo()
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
        solversViewer, "%d\t%e\n", nIters, res); CHKERRQ(ierr);

    ierr = PetscLogStagePop(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // writeLinSolversInfo

// write the time value into a HDF5 file
PetscErrorCode NavierStokesSolver::writeTimeHDF5(const PetscReal &t,
                                                 const std::string &filePath)
{
    PetscErrorCode ierr;
    PetscViewer viewer;

    PetscFunctionBeginUser;

    ierr = PetscViewerCreate(comm, &viewer); CHKERRQ(ierr);
    ierr = PetscViewerSetType(viewer, PETSCVIEWERHDF5); CHKERRQ(ierr);
    ierr = PetscViewerFileSetMode(viewer, FILE_MODE_APPEND); CHKERRQ(ierr);
    ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);
    // attribute has to belong to an existing dataset (choosing p)
    ierr = PetscViewerHDF5WriteAttribute(
        viewer, "/p", "time", PETSC_DOUBLE, &t); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // writeTimeHDF5

// read the time value from a HDF5 file
PetscErrorCode NavierStokesSolver::readTimeHDF5(const std::string &filePath,
                                                PetscReal &t)
{
    PetscErrorCode ierr;
    PetscViewer viewer;

    PetscFunctionBeginUser;

    ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr);
    ierr = PetscViewerSetType(viewer, PETSCVIEWERHDF5); CHKERRQ(ierr);
    ierr = PetscViewerFileSetMode(viewer, FILE_MODE_READ); CHKERRQ(ierr);
    ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);
    // attribute has to belong to an existing dataset (choosing p)
    ierr = PetscViewerHDF5ReadAttribute(viewer, "/p", "time", PETSC_DOUBLE, &t);
    CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // readTimeHDF5

// monitor the solution at probes
PetscErrorCode NavierStokesSolver::monitorProbes()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = PetscLogStagePush(stageMonitor); CHKERRQ(ierr);

    for (auto probe : probes)
    {
        if (ite % probe->nsave == 0 && t >= probe->tstart && t <= probe->tend)
        {
            ierr = probe->monitor(solution, mesh, t); CHKERRQ(ierr);
        }
    }

    ierr = PetscLogStagePop(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // monitorProbes

