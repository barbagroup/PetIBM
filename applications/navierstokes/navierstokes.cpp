/**
 * \file navierstokes.cpp
 * \brief Implementation of the class \c NavierStokesSolver.
 * \see nssolver
 * \ingroup nssolver
 */

// PETSc
# include <petscviewerhdf5.h>

// Navier-Stokes solver
# include "navierstokes.h"

using namespace petibm;


NavierStokesSolver::NavierStokesSolver(const petibm::type::Mesh &inMesh,
            const petibm::type::Boundary &inBC,const YAML::Node &node)
{
    initialize(inMesh, inBC, node);
} // NavierStokesSolver


NavierStokesSolver::~NavierStokesSolver()
{
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    PetscBool finalized;

    ierr = PetscFinalized(&finalized); CHKERRV(ierr);
    if (finalized) return;

    ierr = VecDestroy(&dP); CHKERRV(ierr);
    ierr = VecDestroy(&bc1); CHKERRV(ierr);
    ierr = VecDestroy(&rhs1); CHKERRV(ierr);
    ierr = VecDestroy(&rhs2); CHKERRV(ierr);
    for (unsigned int i=0; i<conv.size(); ++i) {
        ierr = VecDestroy(&conv[i]); CHKERRV(ierr);
    }
    for (unsigned int i=0; i<diff.size(); ++i) {
        ierr = VecDestroy(&diff[i]); CHKERRV(ierr);
    }

    ierr = MatDestroy(&A); CHKERRV(ierr);
    ierr = MatDestroy(&DBNG); CHKERRV(ierr);
    ierr = MatDestroy(&BNG); CHKERRV(ierr);
    ierr = MatDestroy(&N); CHKERRV(ierr);
    ierr = MatDestroy(&G); CHKERRV(ierr);
    ierr = MatDestroy(&D); CHKERRV(ierr);
    ierr = MatDestroy(&DCorrection); CHKERRV(ierr);
    ierr = MatDestroy(&L); CHKERRV(ierr);
    ierr = MatDestroy(&LCorrection); CHKERRV(ierr);

    for(auto &it: asciiViewers) {
        ierr = PetscViewerDestroy(&it.second); CHKERRV(ierr);
    }
} // ~NavierStokesSolver


PetscErrorCode NavierStokesSolver::destroy()
{
    PetscFunctionBeginUser;
    PetscErrorCode ierr;

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
    
    settings.reset();
    bc.reset(); // decrease reference count or destroy
    solution.reset(); // decrease reference count or destroy
    mesh.reset(); // decrease reference count or destroy
    convCoeffs.reset(); // decrease reference count or destroy
    diffCoeffs.reset(); // decrease reference count or destroy
    vSolver.reset(); // decrease reference count or destroy
    pSolver.reset(); // decrease reference count or destroy

    isRefP = PETSC_FALSE;
    dt = 0.0;
    nu = 0.0;

    for(auto &it: asciiViewers) {
        ierr = PetscViewerDestroy(&it.second); CHKERRQ(ierr);
    }
    asciiViewers.clear();

    PetscFunctionReturn(0);
} // destroy


PetscErrorCode NavierStokesSolver::initialize(const petibm::type::Mesh &inMesh,
            const petibm::type::Boundary &inBC,const YAML::Node &node)
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

    // create PETSc Vecs
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

    ierr = PetscLogStagePop(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // initialize


// create operators, i.e., PETSc Mats
PetscErrorCode NavierStokesSolver::createOperators()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    
    Mat     BN; // a temporary operator
    
    // create basic operators
    ierr = operators::createDivergence(
            mesh, bc, D, DCorrection, PETSC_FALSE); CHKERRQ(ierr);
    ierr = operators::createGradient(mesh, G, PETSC_FALSE); CHKERRQ(ierr);
    ierr = operators::createLaplacian(mesh, bc, L, LCorrection); CHKERRQ(ierr);
    ierr = operators::createConvection(mesh, bc, N); CHKERRQ(ierr);
    
    // create combined operator: A
    ierr = MatDuplicate(L, MAT_COPY_VALUES, &A); CHKERRQ(ierr);
    ierr = MatScale(A, -diffCoeffs->implicitCoeff * nu); CHKERRQ(ierr);
    ierr = MatShift(A, 1.0/dt); CHKERRQ(ierr);
    
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

    PetscFunctionReturn(0);
} // createOperators


PetscErrorCode NavierStokesSolver::createVectors()
{
    PetscFunctionBeginUser;
    
    PetscErrorCode ierr;
    
    ierr = VecDuplicate(solution->pGlobal, &dP); CHKERRQ(ierr);
    ierr = VecDuplicate(solution->UGlobal, &bc1); CHKERRQ(ierr);
    ierr = VecDuplicate(solution->UGlobal, &rhs1); CHKERRQ(ierr);
    ierr = VecDuplicate(solution->pGlobal, &rhs2); CHKERRQ(ierr);
    
    conv.resize(convCoeffs->nExplicit);
    for (unsigned int i=0; i<conv.size(); ++i) {
        ierr = VecDuplicate(solution->UGlobal, &conv[i]); CHKERRQ(ierr);
    }
    
    diff.resize(diffCoeffs->nExplicit);
    for (unsigned int i=0; i<diff.size(); ++i) {
        ierr = VecDuplicate(solution->UGlobal, &diff[i]); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
} // createVectors


PetscErrorCode NavierStokesSolver::setNullSpace()
{
    PetscFunctionBeginUser;
    
    PetscErrorCode  ierr;
    
    std::string type;
    
    ierr = pSolver->getType(type); CHKERRQ(ierr);
    
    if (type == "PETSc KSP")
    {
        MatNullSpace nsp;
        ierr = MatNullSpaceCreate(
                mesh->comm, PETSC_TRUE, 0, nullptr, &nsp); CHKERRQ(ierr);
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
} // setNullSpace


// advance the flow solver for one time-step
PetscErrorCode NavierStokesSolver::advance()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    // prepare velocity system and solve it
    ierr = assembleRHSVelocity(); CHKERRQ(ierr);
    ierr = solveVelocity(); CHKERRQ(ierr);

    // prepare Poisson system and solve it
    ierr = assembleRHSPoisson(); CHKERRQ(ierr);
    ierr = solvePoisson(); CHKERRQ(ierr);

    // correct solutions
    ierr = projectionStep(); CHKERRQ(ierr);

    // update values of ghost points
    ierr = bc->updateGhostValues(solution); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // advance


// prepare right-hand-side vector for velocity solver
PetscErrorCode NavierStokesSolver::assembleRHSVelocity()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = PetscLogStagePush(stageRHSVelocity); CHKERRQ(ierr);

    // initialize RHS with pressure gradient at time-step n
    ierr = MatMult(G, solution->pGlobal, rhs1); CHKERRQ(ierr);
    ierr = VecScale(rhs1, -1.0); CHKERRQ(ierr);

    // add explicit part of time derivative ($\frac{u^n}{\Delta t}$)
    ierr = VecAXPY(rhs1, 1.0/dt, solution->UGlobal); CHKERRQ(ierr);

    // add convection terms from time index n, n-1, n-2, ... to rhs1:
    {
        // 1. discard the oldest-time-step term, and decrease the time-step by 1
        for(int i=conv.size()-1; i>0; i--) {
            ierr = VecSwap(conv[i], conv[i-1]); CHKERRQ(ierr);
        }
        
        if (conv.size() > 0)
        {
            // 2. get convection term at time-step n
            ierr = MatMult(N, solution->UGlobal, conv[0]); CHKERRQ(ierr);
        
            // 3. move n-th time-step convection term to right-hand-side
            ierr = VecScale(conv[0], -1.0); CHKERRQ(ierr);
        }
        
        // 4. add all explicit convective terms to rhs1
        for(unsigned int i=0; i<conv.size(); ++i) {
            ierr = VecAXPY(rhs1, convCoeffs->explicitCoeffs[i], conv[i]);
            CHKERRQ(ierr);
        }
    }

    // add diffusion terms from time index n, n-1, n-2, ... to rhs1:
    {
        // 1. discard the oldest-time-step term, and decrease the time-step by 1
        for(int i=diff.size()-1; i>0; i--) {
            ierr = VecSwap(diff[i], diff[i-1]); CHKERRQ(ierr);
        }
        
        // 2. get diffusion term at time-step n
        if (diff.size() > 0)
        {
            ierr = MatMult(L, solution->UGlobal, diff[0]); CHKERRQ(ierr);
            ierr = MatMultAdd(LCorrection, 
                    solution->UGlobal, diff[0], diff[0]); CHKERRQ(ierr);
            ierr = VecScale(diff[0], nu); CHKERRQ(ierr);
        }
        
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
PetscErrorCode NavierStokesSolver::solveVelocity()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = PetscLogStagePush(stageSolveVelocity); CHKERRQ(ierr);

    ierr = vSolver->solve(solution->UGlobal, rhs1); CHKERRQ(ierr);

    ierr = PetscLogStagePop(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // solveVelocity


// prepare tight-hand-side vector for Poisson system
PetscErrorCode NavierStokesSolver::assembleRHSPoisson()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = PetscLogStagePush(stageRHSPoisson); CHKERRQ(ierr);

    // get Du* (which is equal to (D_{interior} + D_{boundary})u*
    ierr = MatMult(D, solution->UGlobal, rhs2); CHKERRQ(ierr);
    ierr = MatMultAdd(DCorrection, solution->UGlobal, rhs2, rhs2); CHKERRQ(ierr);
    
    
    if (isRefP)
    {
        ierr = VecSetValue(rhs2, 0, 0.0, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecAssemblyBegin(rhs2); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(rhs2); CHKERRQ(ierr);
    }

    ierr = PetscLogStagePop(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // assembleRHSPoisson


// solve Poisson system
PetscErrorCode NavierStokesSolver::solvePoisson()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = PetscLogStagePush(stageSolvePoisson); CHKERRQ(ierr);

    // solve for the pressure correction
    ierr = pSolver->solve(dP, rhs2); CHKERRQ(ierr);

    ierr = PetscLogStagePop(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // solvePoisson


// projection
PetscErrorCode NavierStokesSolver::projectionStep()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = PetscLogStagePush(stageProjectionStep); CHKERRQ(ierr);

    // project velocity field onto divergence-free space
    ierr = MatMult(BNG, dP, rhs1); CHKERRQ(ierr);
    ierr = VecAXPY(solution->UGlobal, -1.0, rhs1); CHKERRQ(ierr);
    
    // add pressure increment
    ierr = VecAXPY(solution->pGlobal, 1.0, dP); CHKERRQ(ierr);

    ierr = PetscLogStagePop(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // projectionStep


// output solutions to the user provided file
PetscErrorCode NavierStokesSolver::write(const std::string &filePath)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = PetscLogStagePush(stageWrite); CHKERRQ(ierr);

    ierr = solution->write(filePath); CHKERRQ(ierr);

    ierr = PetscLogStagePop(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // write


// output extra data required for restarting to the user provided file
PetscErrorCode NavierStokesSolver::writeRestartData(const std::string &filePath)
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
} // writeRestartData


// read data necessary for restarting
PetscErrorCode NavierStokesSolver::readRestartData(const std::string &filePath)
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
} // readRestartData


PetscErrorCode NavierStokesSolver::initializeASCIIFiles(
        const std::string &filePath, const PetscFileMode &mode)
{
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    // create ASCII viewer
    ierr = PetscViewerCreate(PETSC_COMM_WORLD, &asciiViewers[filePath]); CHKERRQ(ierr);
    ierr = PetscViewerSetType(asciiViewers[filePath], PETSCVIEWERASCII); CHKERRQ(ierr);
    ierr = PetscViewerFileSetMode(asciiViewers[filePath], mode); CHKERRQ(ierr);
    ierr = PetscViewerFileSetName(asciiViewers[filePath], filePath.c_str()); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // initializeASCIIFiles


// write numbers of iterations and residuals of linear solvers to a file
PetscErrorCode NavierStokesSolver::writeIterations(
        const int &timeIndex, const std::string &filePath)
{
    PetscErrorCode ierr;
    PetscInt nIters;
    PetscReal res;

    PetscFunctionBeginUser;
    
    // write current time
    ierr = PetscViewerASCIIPrintf(
            asciiViewers[filePath], "%d", timeIndex); CHKERRQ(ierr);
    
    ierr = vSolver->getIters(nIters); CHKERRQ(ierr);
    ierr = vSolver->getResidual(res); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(
            asciiViewers[filePath], "\t%d\t%e", nIters, res); CHKERRQ(ierr);
    
    ierr = pSolver->getIters(nIters); CHKERRQ(ierr);
    ierr = pSolver->getResidual(res); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(
            asciiViewers[filePath], "\t%d\t%e\n", nIters, res); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
} // writeIterations
