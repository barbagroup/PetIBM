/**
 * \file navierstokes.cpp
 * \brief Implementation of the class \c TairaColoniusSolver.
 */

// PETSc
# include <petscviewerhdf5.h>

// Navier-Stokes solver
# include "tairacolonius.h"

using namespace petibm;


TairaColoniusSolver::TairaColoniusSolver(
        const petibm::type::Mesh &inMesh, const petibm::type::Boundary &inBC,
        const petibm::type::BodyPack &inBodies,const YAML::Node &node)
{
    initialize(inMesh, inBC, inBodies, node);
} // TairaColoniusSolver


PetscErrorCode TairaColoniusSolver::initialize(
        const petibm::type::Mesh &inMesh, const petibm::type::Boundary &inBC,
        const petibm::type::BodyPack &inBodies,const YAML::Node &node)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = PetscLogStageRegister("initialize", &stageInitialize); CHKERRQ(ierr);
    ierr = PetscLogStagePush(stageInitialize); CHKERRQ(ierr);
    
    // save the reference to the YAML::Node
    settings = node;
    
    // save the reference to mesh and boundary conditions
    mesh = inMesh;
    bc = inBC;
    bodies = inBodies;
    
    // get dt and nu
    // TODO: should we check the keys "parameters", "dt", "flow", and "nu"?
    dt = settings["parameters"]["dt"].as<PetscReal>();
    nu = settings["flow"]["nu"].as<PetscReal>();
    
    // create solution object
    ierr = solution::createSolution(mesh, solution); CHKERRQ(ierr);

    // create time schemes objects
    ierr = timeintegration::createTimeIntegration("convection", settings, convCoeffs); CHKERRQ(ierr);
    ierr = timeintegration::createTimeIntegration("diffusion", settings, diffCoeffs); CHKERRQ(ierr);
    
    // create linear solve objects
    ierr = linsolver::createLinSolver("velocity", settings, vSolver); CHKERRQ(ierr);
    ierr = linsolver::createLinSolver("poisson", settings, pSolver); CHKERRQ(ierr);
    
    // set initial values to solutions
    ierr = solution->applyIC(settings); CHKERRQ(ierr);

    // initialize ghost points values and Eqs; must before creating operators
    ierr = bc->setGhostICs(solution); CHKERRQ(ierr);

    // create operators (PETSc Mats)
    ierr = createOperators(); CHKERRQ(ierr);

    // creater PETSc Vecs
    ierr = createVectors(); CHKERRQ(ierr);

    // set coefficient matrix to linear solvers
    ierr = vSolver->setMatrix(A); CHKERRQ(ierr);
    ierr = pSolver->setMatrix(DBNG); CHKERRQ(ierr);
    
    // register events
    ierr = PetscLogStageRegister("rhsVelocity", &stageRHSVelocity); CHKERRQ(ierr);
    ierr = PetscLogStageRegister("solveVelocity", &stageSolveVelocity); CHKERRQ(ierr);
    ierr = PetscLogStageRegister("rhsPoisson", &stageRHSPoisson); CHKERRQ(ierr);
    ierr = PetscLogStageRegister("solvePoisson", &stageSolvePoisson); CHKERRQ(ierr);
    ierr = PetscLogStageRegister("projectionStep", &stageProjectionStep); CHKERRQ(ierr);
    ierr = PetscLogStageRegister("write", &stageWrite); CHKERRQ(ierr);
    ierr = PetscLogStageRegister("integrateForces", &stageIntegrateForces); CHKERRQ(ierr);

    ierr = PetscLogStagePop(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // initialize


// create operators, i.e., PETSc Mats
PetscErrorCode TairaColoniusSolver::createOperators()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    
    Mat     R, MHat, BN, GH[2], DE[2]; // a temporary operators
    Vec     RDiag, MHatDiag; // a temporary vectors
    IS      is;
    
    // create basic operators
    ierr = operators::createDivergence(
            mesh, bc, DE[0], DCorrection, PETSC_FALSE); CHKERRQ(ierr);
    ierr = operators::createGradient(mesh, GH[0], PETSC_FALSE); CHKERRQ(ierr);
    ierr = operators::createLaplacian(mesh, bc, L, LCorrection); CHKERRQ(ierr);
    ierr = operators::createConvection(mesh, bc, N); CHKERRQ(ierr);
    
    // create combined operator: A
    ierr = MatDuplicate(L, MAT_COPY_VALUES, &A); CHKERRQ(ierr);
    ierr = MatScale(A, -diffCoeffs->implicitCoeff * nu); CHKERRQ(ierr);
    ierr = MatShift(A, 1.0/dt); CHKERRQ(ierr);

    // create diagonal matrix R and get a Vec holding the diagonal
    ierr = petibm::operators::createR(mesh, R); CHKERRQ(ierr);
    ierr = MatCreateVecs(R, nullptr, &RDiag); CHKERRQ(ierr);
    ierr = MatGetDiagonal(R, RDiag); CHKERRQ(ierr);
    
    // create diagonal matrix M and get a Vec holding the diagonal
    ierr = petibm::operators::createMHead(mesh, MHat); CHKERRQ(ierr);
    ierr = MatCreateVecs(MHat, nullptr, &MHatDiag); CHKERRQ(ierr);
    ierr = MatGetDiagonal(MHat, MHatDiag); CHKERRQ(ierr);
    
    // create a Delta operator and its transpose (equal to H)
    ierr = petibm::operators::createDelta(mesh, bc, bodies, DE[1]); CHKERRQ(ierr);
    ierr = MatTranspose(DE[1], MAT_INITIAL_MATRIX, &GH[1]); CHKERRQ(ierr);
    
    // create operator E
    ierr = MatDiagonalScale(DE[1], nullptr, RDiag); CHKERRQ(ierr);
    ierr = MatDiagonalScale(DE[1], nullptr, MHatDiag); CHKERRQ(ierr);
    ierr = VecDestroy(&RDiag); CHKERRQ(ierr);
    ierr = VecDestroy(&MHatDiag); CHKERRQ(ierr);
    ierr = MatDestroy(&R); CHKERRQ(ierr);
    ierr = MatDestroy(&MHat); CHKERRQ(ierr);
    
    // create -H
    ierr = MatScale(GH[1], -1.0); CHKERRQ(ierr);

    // get combined operator; R is used temporarily
    ierr = MatCreateNest(mesh->comm, 1, nullptr, 2, nullptr, GH, &R); CHKERRQ(ierr);
    ierr = MatConvert(R, MATAIJ, MAT_INITIAL_MATRIX, &G); CHKERRQ(ierr);
    ierr = MatDestroy(&R); CHKERRQ(ierr);
    
    IS iss[2];
    ierr = MatCreateNest(mesh->comm, 2, nullptr, 1, nullptr, DE, &R); CHKERRQ(ierr);
    ierr = MatConvert(R, MATAIJ, MAT_INITIAL_MATRIX, &D); CHKERRQ(ierr);
    ierr = MatNestGetISs(R, iss, nullptr); CHKERRQ(ierr);
    ierr = ISDuplicate(iss[0], &isDE[0]); CHKERRQ(ierr);
    ierr = ISDuplicate(iss[1], &isDE[1]); CHKERRQ(ierr);
    ierr = ISCopy(iss[0], isDE[0]); CHKERRQ(ierr);
    ierr = ISCopy(iss[1], isDE[1]); CHKERRQ(ierr);
    ierr = MatDestroy(&R); CHKERRQ(ierr);
    
    // create combined operator: BNG
    ierr = operators::createBnHead(
            L, dt, diffCoeffs->implicitCoeff * nu, 1, BN); CHKERRQ(ierr);
    ierr = MatMatMult(
            BN, G, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &BNG); CHKERRQ(ierr);
    
    // create combined operator: DBNG
    ierr = MatMatMult(
            D, BNG, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &DBNG); CHKERRQ(ierr);
    
    // set null space to DBNG
    ierr = setNullSpace(); CHKERRQ(ierr);
    
    // destroy temporary operator
    ierr = MatDestroy(&BN); CHKERRQ(ierr);
    
    // register auto-destroy objects
    ierr = PetscObjectRegisterDestroy((PetscObject) DBNG); CHKERRQ(ierr);
    ierr = PetscObjectRegisterDestroy((PetscObject) BNG); CHKERRQ(ierr);
    ierr = PetscObjectRegisterDestroy((PetscObject) A); CHKERRQ(ierr);
    ierr = PetscObjectRegisterDestroy((PetscObject) N); CHKERRQ(ierr);
    ierr = PetscObjectRegisterDestroy((PetscObject) LCorrection); CHKERRQ(ierr);
    ierr = PetscObjectRegisterDestroy((PetscObject) L); CHKERRQ(ierr);
    ierr = PetscObjectRegisterDestroy((PetscObject) G); CHKERRQ(ierr);
    ierr = PetscObjectRegisterDestroy((PetscObject) DCorrection); CHKERRQ(ierr);
    ierr = PetscObjectRegisterDestroy((PetscObject) D); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // assembleOperators


PetscErrorCode TairaColoniusSolver::createVectors()
{
    PetscFunctionBeginUser;
    
    PetscErrorCode ierr;
    
    ierr = VecDuplicate(solution->UGlobal, &bc1); CHKERRQ(ierr);
    ierr = VecDuplicate(solution->UGlobal, &rhs1); CHKERRQ(ierr);
    
    conv.resize(convCoeffs->nExplicit);
    for (unsigned int i=0; i<conv.size(); ++i) {
        ierr = VecDuplicate(solution->UGlobal, &conv[i]); CHKERRQ(ierr);
    }
    
    diff.resize(diffCoeffs->nExplicit);
    for (unsigned int i=0; i<diff.size(); ++i) {
        ierr = VecDuplicate(solution->UGlobal, &diff[i]); CHKERRQ(ierr);
    }
    
    // phi, dphi, rhs2
    ierr = MatCreateVecs(DBNG, &dphi, &rhs2); CHKERRQ(ierr);
    ierr = VecDuplicate(dphi, &phi); CHKERRQ(ierr);
    
    // reset pGlobal
    //ierr = VecDestroy(&(solution->pGlobal)); CHKERRQ(ierr);
    //ierr = VecGetSubVector(phi, isDE[0], &(solution->pGlobal)); CHKERRQ(ierr);
    Vec temp;
    PetscReal *a;
    ierr = VecGetSubVector(phi, isDE[0], &temp); CHKERRQ(ierr);
    ierr = VecGetArray(temp, &a); CHKERRQ(ierr);
    ierr = VecPlaceArray(solution->pGlobal, a); CHKERRQ(ierr);
    a = nullptr;
    
    // f
    ierr = VecGetSubVector(phi, isDE[1], &f); CHKERRQ(ierr);
    
    // bc2
    ierr = VecGetSubVector(rhs2, isDE[0], &bc2); CHKERRQ(ierr);
    
    // resister auto-destroy
    ierr = PetscObjectRegisterDestroy((PetscObject) rhs1); CHKERRQ(ierr);
    ierr = PetscObjectRegisterDestroy((PetscObject) bc1); CHKERRQ(ierr);
    ierr = PetscObjectRegisterDestroy((PetscObject) phi); CHKERRQ(ierr);
    ierr = PetscObjectRegisterDestroy((PetscObject) dphi); CHKERRQ(ierr);
    
    for (unsigned int i=0; i<conv.size(); ++i) {
        ierr = PetscObjectRegisterDestroy((PetscObject) conv[i]); CHKERRQ(ierr);
    }
    
    for (unsigned int i=0; i<diff.size(); ++i) {
        ierr = PetscObjectRegisterDestroy((PetscObject) diff[i]); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}


PetscErrorCode TairaColoniusSolver::setNullSpace()
{
    PetscFunctionBeginUser;
    
    PetscErrorCode  ierr;
    
    std::string type;
    
    ierr = pSolver->getType(type); CHKERRQ(ierr);
    
    if (type == "PETSc KSP")
    {
        Vec n, phiPortion;
        ierr = MatCreateVecs(DBNG, &n, nullptr); CHKERRQ(ierr);
        ierr = VecSet(n, 0.0); CHKERRQ(ierr);
        ierr = VecGetSubVector(n, isDE[0], &phiPortion); CHKERRQ(ierr);
        ierr = VecSet(phiPortion, 1.0 / std::sqrt(mesh->pN)); CHKERRQ(ierr);
        ierr = VecRestoreSubVector(n, isDE[0], &phiPortion); CHKERRQ(ierr);
        MatNullSpace nsp;
        ierr = MatNullSpaceCreate(
                mesh->comm, PETSC_FALSE, 1, &n, &nsp); CHKERRQ(ierr);
        ierr = MatSetNullSpace(DBNG, nsp); CHKERRQ(ierr);
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
}


// advance the flow solver for one time-step
PetscErrorCode TairaColoniusSolver::advance()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    // prepare velocity system and solve it
    ierr = assembleRHSVelocity(); CHKERRQ(ierr);
    ierr = solveVelocity(); CHKERRQ(ierr);

    // prepare poisson system and solve it
    ierr = assembleRHSPoisson(); CHKERRQ(ierr);
    ierr = solvePoisson(); CHKERRQ(ierr);

    // correct solutions
    ierr = projectionStep(); CHKERRQ(ierr);

    // update values of ghost points
    ierr = bc->updateGhostValues(solution); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // solve


// prepare right-hand-side vector for velocity solver
PetscErrorCode TairaColoniusSolver::assembleRHSVelocity()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = PetscLogStagePush(stageRHSVelocity); CHKERRQ(ierr);

    // initialize RHS with pressure gradient at time-step n
    ierr = MatMult(G, phi, rhs1); CHKERRQ(ierr);
    ierr = VecScale(rhs1, -1.0); CHKERRQ(ierr);

    // add explicit part of time derivative ($\frac{u^n}{\Delta t}$)
    ierr = VecAXPY(rhs1, 1.0/dt, solution->UGlobal); CHKERRQ(ierr);

    // add convection terms from time index n, n-1, n-2, ... to rhs1:
    {
        // 1. discard the oldest-time-step term, and decrease the time-step by 1
        for(unsigned int i=conv.size()-1; i>0; i--) {
            ierr = VecSwap(conv[i], conv[i-1]); CHKERRQ(ierr);
        }
        
        // 2. get convection term at time-step n
        ierr = MatMult(N, solution->UGlobal, conv[0]); CHKERRQ(ierr);
        
        // 3. move n-th time-step convection term to right-hand-side
        ierr = VecScale(conv[0], -1.0); CHKERRQ(ierr);
        
        // 4. add all explicit convective terms to rhs1
        for(unsigned int i=0; i<conv.size(); ++i) {
            ierr = VecAXPY(rhs1, convCoeffs->explicitCoeffs[i], conv[i]);
            CHKERRQ(ierr);
        }
    }

    // add diffusion terms from time index n, n-1, n-2, ... to rhs1:
    {
        // 1. discard the oldest-time-step term, and decrease the time-step by 1
        for(unsigned int i=diff.size()-1; i>0; i--) {
            ierr = VecSwap(diff[i], diff[i-1]); CHKERRQ(ierr);
        }
        
        // 2. get diffusion term at time-step n
        ierr = MatMult(L, solution->UGlobal, diff[0]); CHKERRQ(ierr);
        ierr = MatMultAdd(LCorrection, 
                solution->UGlobal, diff[0], diff[0]); CHKERRQ(ierr);
        ierr = VecScale(diff[0], nu); CHKERRQ(ierr);
        
        // 3. add all explicit diffusive terms to rhs1
        for(unsigned int i=0; i<diff.size(); ++i) {
            ierr = VecAXPY(rhs1, diffCoeffs->explicitCoeffs[i], diff[i]);
            CHKERRQ(ierr);
        }
    }
    
    // add the diffusion BC correction from time index n+1 to rhs1
    {
        // 1. update the Eq.s of ghost points to time index n+1
        ierr = bc->updateEqs(solution, dt); CHKERRQ(ierr);
        
        // 2. get BC correction of Laplacian of time n+1
        ierr = MatMult(LCorrection, solution->UGlobal, bc1); CHKERRQ(ierr);
        ierr = VecScale(bc1, nu); CHKERRQ(ierr);
        
        // 3. add bc1 to rhs
        ierr = VecAXPY(rhs1, diffCoeffs->implicitCoeff, bc1); CHKERRQ(ierr);
    }
    
    ierr = PetscLogStagePop(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // assembleRHSVelocity


// solve velocity system
PetscErrorCode TairaColoniusSolver::solveVelocity()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = PetscLogStagePush(stageSolveVelocity); CHKERRQ(ierr);

    ierr = vSolver->solve(solution->UGlobal, rhs1); CHKERRQ(ierr);

    ierr = PetscLogStagePop(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // solveVelocity


// prepare tight-hand-side vector for Poisson system
PetscErrorCode TairaColoniusSolver::assembleRHSPoisson()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = PetscLogStagePush(stageRHSPoisson); CHKERRQ(ierr);

    // get Du* (which is equal to (D_{interior} + D_{boundary})u*
    // note, bc2 is a subset of rhs2
    ierr = VecSet(rhs2, 0.0); CHKERRQ(ierr);
    ierr = MatMult(DCorrection, solution->UGlobal, bc2); CHKERRQ(ierr);
    ierr = MatMultAdd(D, solution->UGlobal, rhs2, rhs2); CHKERRQ(ierr);
    
    
    if (isRefP)
    {
        ierr = VecSetValue(rhs2, 0, 0.0, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecAssemblyBegin(rhs2); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(rhs2); CHKERRQ(ierr);
    }
    
    ierr = PetscLogStagePop(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // assembleRHSPoisson


// solve poisson system
PetscErrorCode TairaColoniusSolver::solvePoisson()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = PetscLogStagePush(stageSolvePoisson); CHKERRQ(ierr);

    // solve for the pressure correction
    ierr = pSolver->solve(dphi, rhs2); CHKERRQ(ierr);

    ierr = PetscLogStagePop(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // solvePoisson


// projection
PetscErrorCode TairaColoniusSolver::projectionStep()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = PetscLogStagePush(stageProjectionStep); CHKERRQ(ierr);

    // project velocity field onto divergence-free space
    ierr = MatMult(BNG, dphi, rhs1); CHKERRQ(ierr);
    ierr = VecAXPY(solution->UGlobal, -1.0, rhs1); CHKERRQ(ierr);
    
    // add phi increment; pressure and forces are subsets (i.e., included)
    ierr = VecAXPY(phi, 1.0, dphi); CHKERRQ(ierr);

    ierr = PetscLogStagePop(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // projectionStep


// output solutions to the user provided file
PetscErrorCode TairaColoniusSolver::write(const std::string &filePath)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = PetscLogStagePush(stageWrite); CHKERRQ(ierr);

    ierr = solution->write(filePath); CHKERRQ(ierr);

    ierr = PetscLogStagePop(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // write


// output extra data required for restarting to the user provided file
PetscErrorCode TairaColoniusSolver::writeRestartData(const std::string &filePath)
{
    PetscFunctionBeginUser;
    
    PetscErrorCode ierr;
    
    PetscBool   fileExist = PETSC_FALSE;
    
    PetscViewer     viewer;

    ierr = PetscLogStagePush(stageWrite); CHKERRQ(ierr);
    
    // test if file exist 
    ierr = PetscTestFile((filePath+".h5").c_str(), 'w', &fileExist); CHKERRQ(ierr);
    
    if (! fileExist) // if not, create one and write u, v, w, and p into it
    {
        ierr = solution->write(filePath); CHKERRQ(ierr);
    }
    // TODO: should we check if the file exist but data is not up-to-date?
    
    // create PetscViewer with append mode
    ierr = PetscViewerCreate(mesh->comm, &viewer); CHKERRQ(ierr);
    ierr = PetscViewerSetType(viewer, PETSCVIEWERHDF5); CHKERRQ(ierr);
    ierr = PetscViewerFileSetMode(viewer, FILE_MODE_APPEND); CHKERRQ(ierr);
    ierr = PetscViewerFileSetName(viewer, (filePath+".h5").c_str()); CHKERRQ(ierr);
    
    // go to the root node first. Not necessary. Just in case.
    ierr = PetscViewerHDF5PushGroup(viewer, "/"); CHKERRQ(ierr);
    
    // write extra data to the end of the file
    // save explicit convection terms first
    ierr = PetscViewerHDF5PushGroup(viewer, "/convection"); CHKERRQ(ierr);
    for(unsigned int i=0; i<conv.size(); ++i)
    {
        ierr = PetscObjectSetName((PetscObject) conv[i],
                std::to_string(i).c_str()); CHKERRQ(ierr);
        ierr = VecView(conv[i], viewer); CHKERRQ(ierr);
    }
    
    // then save explicit diffusion terms
    ierr = PetscViewerHDF5PushGroup(viewer, "/diffusion"); CHKERRQ(ierr);
    for(unsigned int i=0; i<diff.size(); ++i)
    {
        ierr = PetscObjectSetName((PetscObject) diff[i],
                std::to_string(i).c_str()); CHKERRQ(ierr);
        ierr = VecView(diff[i], viewer); CHKERRQ(ierr);
    }
    
    // destroy viewer
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

    ierr = PetscLogStagePop(); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}


// read data necessary for restarting
PetscErrorCode TairaColoniusSolver::readRestartData(const std::string &filePath)
{
    PetscFunctionBeginUser;
    
    PetscErrorCode  ierr;
    
    PetscBool   fileExist = PETSC_FALSE;
    
    PetscViewer     viewer;

    // test if file exist 
    ierr = PetscTestFile((filePath+".h5").c_str(), 'r', &fileExist); CHKERRQ(ierr);
    
    if (! fileExist) // if not, return error
        SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_FILE_READ,
                "Could not find file \"%s\" for restarting.\n", 
                (filePath + ".h5").c_str());
    
    // read primary fields
    ierr = solution->read(filePath); CHKERRQ(ierr);
    
    // create PetscViewer with append mode
    ierr = PetscViewerCreate(mesh->comm, &viewer); CHKERRQ(ierr);
    ierr = PetscViewerSetType(viewer, PETSCVIEWERHDF5); CHKERRQ(ierr);
    ierr = PetscViewerFileSetMode(viewer, FILE_MODE_READ); CHKERRQ(ierr);
    ierr = PetscViewerFileSetName(viewer, (filePath+".h5").c_str()); CHKERRQ(ierr);
    
    // go to the root node first. Not necessary. Just in case.
    ierr = PetscViewerHDF5PushGroup(viewer, "/"); CHKERRQ(ierr);
    
    // read explicit convection terms first
    ierr = PetscViewerHDF5PushGroup(viewer, "/convection"); CHKERRQ(ierr);
    for(unsigned int i=0; i<conv.size(); ++i)
    {
        ierr = PetscObjectSetName((PetscObject) conv[i],
                std::to_string(i).c_str()); CHKERRQ(ierr);
        ierr = VecLoad(conv[i], viewer); CHKERRQ(ierr);
    }
    
    // then read explicit diffusion terms
    ierr = PetscViewerHDF5PushGroup(viewer, "/diffusion"); CHKERRQ(ierr);
    for(unsigned int i=0; i<diff.size(); ++i)
    {
        ierr = PetscObjectSetName((PetscObject) diff[i],
                std::to_string(i).c_str()); CHKERRQ(ierr);
        ierr = VecLoad(diff[i], viewer); CHKERRQ(ierr);
    }
    
    // destroy viewer
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
    
    // update values and Eqs of ghost points based on current solution
    // TODO: for convective BCs, it's not totally correct
    ierr = bc->setGhostICs(solution); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}


// write numbers of iterations and residuals of linear solvers to a file
PetscErrorCode TairaColoniusSolver::writeIterations(
        const int &timeIndex, const std::string &filePath)
{
    PetscErrorCode ierr;
    PetscInt nIters;
    PetscReal res;
    PetscViewer viewer;
    static PetscFileMode mode = FILE_MODE_WRITE;

    PetscFunctionBeginUser;
    
    // create ASCII viewer
    ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr);
    ierr = PetscViewerSetType(viewer, PETSCVIEWERASCII); CHKERRQ(ierr);
    ierr = PetscViewerFileSetMode(viewer, mode); CHKERRQ(ierr);
    ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);
    
    // write current time
    ierr = PetscViewerASCIIPrintf(viewer, "%d", timeIndex); CHKERRQ(ierr);
    
    ierr = vSolver->getIters(nIters); CHKERRQ(ierr);
    ierr = vSolver->getResidual(res); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer, "\t%d\t%e", nIters, res); CHKERRQ(ierr);
    
    ierr = pSolver->getIters(nIters); CHKERRQ(ierr);
    ierr = pSolver->getResidual(res); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer, "\t%d\t%e\n", nIters, res); CHKERRQ(ierr);
    
    // destroy
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
    
    // next time we'll append data
    mode = FILE_MODE_APPEND;

    PetscFunctionReturn(0);
} // writeIterations


// write averaged forces
PetscErrorCode TairaColoniusSolver::writeIntegratedForces(
            const PetscReal &t, const std::string &filePath)
{
    PetscErrorCode ierr;
    PetscViewer         viewer;
    petibm::type::RealVec2D     fAvg;
    static PetscFileMode mode = FILE_MODE_WRITE;

    PetscFunctionBeginUser;

    ierr = PetscLogStagePush(stageIntegrateForces); CHKERRQ(ierr);

    // get averaged forces first
    ierr = bodies->calculateAvgForces(f, fAvg); CHKERRQ(ierr);

    ierr = PetscLogStagePop(); CHKERRQ(ierr);
    
    // create ASCII viewer
    ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr);
    ierr = PetscViewerSetType(viewer, PETSCVIEWERASCII); CHKERRQ(ierr);
    ierr = PetscViewerFileSetMode(viewer, mode); CHKERRQ(ierr);
    ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);
    
    // write current time
    ierr = PetscViewerASCIIPrintf(viewer, "%10.8e", t); CHKERRQ(ierr);
    
    // write forces body by body
    for(unsigned int i=0; i<bodies->nBodies; ++i)
    {
        for(unsigned int d=0; d<mesh->dim; ++d)
        {
            ierr = PetscViewerASCIIPrintf(
                    viewer, "\t%10.8e", fAvg[i][d]); CHKERRQ(ierr);
        }
    }
    ierr = PetscViewerASCIIPrintf(viewer, "\n"); CHKERRQ(ierr);
    
    // destroy
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
    
    // next we'll just append data
    mode = FILE_MODE_APPEND;

    PetscFunctionReturn(0);
} // writeIntegratedForces


// manual finalization
PetscErrorCode TairaColoniusSolver::finalize()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = VecDestroy(&dP); CHKERRQ(ierr);
    ierr = VecDestroy(&bc1); CHKERRQ(ierr);
    ierr = VecDestroy(&rhs1); CHKERRQ(ierr);
    ierr = VecDestroy(&rhs2); CHKERRQ(ierr);
    for (unsigned int i=0; i<conv.size(); ++i) {
        ierr = VecDestroy(&conv[i]); CHKERRQ(ierr);
    }
    for (unsigned int i=0; i<diff.size(); ++i) {
        ierr = VecDestroy(&diff[i]); CHKERRQ(ierr);
    }

    ierr = MatDestroy(&A); CHKERRQ(ierr);
    ierr = MatDestroy(&DBNG); CHKERRQ(ierr);
    ierr = MatDestroy(&BNG); CHKERRQ(ierr);
    ierr = MatDestroy(&N); CHKERRQ(ierr);
    ierr = MatDestroy(&G); CHKERRQ(ierr);
    ierr = MatDestroy(&D); CHKERRQ(ierr);
    ierr = MatDestroy(&DCorrection); CHKERRQ(ierr);
    ierr = MatDestroy(&L); CHKERRQ(ierr);
    ierr = MatDestroy(&LCorrection); CHKERRQ(ierr);
    
    settings.~Node();
    mesh.~shared_ptr(); // decrease reference count or destroy
    bc.~shared_ptr(); // decrease reference count or destroy
    solution.~shared_ptr(); // decrease reference count or destroy
    convCoeffs.~shared_ptr(); // decrease reference count or destroy
    diffCoeffs.~shared_ptr(); // decrease reference count or destroy
    vSolver.~shared_ptr(); // decrease reference count or destroy
    pSolver.~shared_ptr(); // decrease reference count or destroy

    PetscFunctionReturn(0);
} // finalize
