/**
 * \file navierstokes.cpp
 * \brief Implementation of the class \c NavierStokesSolver.
 */

// STL
#include <fstream>
#include <iomanip>
#include <string>

// PETSc
# include <petscviewerhdf5.h>

#include "navierstokes.h"

using namespace petibm;


NavierStokesSolver::NavierStokesSolver(const petibm::type::Mesh &inMesh,
            const petibm::type::Boundary &inBC,const YAML::Node &node)
{
    initialize(inMesh, inBC, node);
} // NavierStokesSolver


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
    ierr = timeintegration::createTimeIntegration("", settings, convCoeffs); CHKERRQ(ierr);
    ierr = timeintegration::createTimeIntegration("", settings, diffCoeffs); CHKERRQ(ierr);
    
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
    PetscLogStageRegister("rhsVelocity", &stageRHSVelocity);
    PetscLogStageRegister("solveVelocity", &stageSolveVelocity);
    PetscLogStageRegister("rhsPoisson", &stageRHSPoisson);
    PetscLogStageRegister("solvePoisson", &stageSolvePoisson);
    PetscLogStageRegister("projectionStep", &stageProjectionStep);
    PetscLogStageRegister("write", &stageWrite);

    ierr = PetscLogStagePop();

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
    {
        MatNullSpace nsp;
        ierr = MatNullSpaceCreate(
                mesh->comm, PETSC_TRUE, 0, nullptr, &nsp); CHKERRQ(ierr);
        ierr = MatSetNullSpace(DBNG, nsp); CHKERRQ(ierr);
        ierr = MatNullSpaceDestroy(&nsp); CHKERRQ(ierr);
    }
    
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
    
    // resister auto-destroy
    ierr = PetscObjectRegisterDestroy((PetscObject) rhs2); CHKERRQ(ierr);
    ierr = PetscObjectRegisterDestroy((PetscObject) rhs1); CHKERRQ(ierr);
    ierr = PetscObjectRegisterDestroy((PetscObject) bc1); CHKERRQ(ierr);
    ierr = PetscObjectRegisterDestroy((PetscObject) dP); CHKERRQ(ierr);
    
    for (unsigned int i=0; i<conv.size(); ++i) {
        ierr = PetscObjectRegisterDestroy((PetscObject) conv[i]); CHKERRQ(ierr);
    }
    
    for (unsigned int i=0; i<conv.size(); ++i) {
        ierr = PetscObjectRegisterDestroy((PetscObject) diff[i]); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}


// advance the flow solver for one time-step
PetscErrorCode NavierStokesSolver::advance()
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

    ierr = PetscLogStagePop(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // assembleRHSPoisson


// solve poisson system
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
    
    // write extra data to the end of the file
    // save explicit convection terms first
    ierr = PetscViewerHDF5PushGroup(viewer, "convection"); CHKERRQ(ierr);
    for(unsigned int i=0; i<convCoeffs->nExplicit; ++i)
    {
        ierr = VecView(conv[i], viewer); CHKERRQ(ierr);
    }
    ierr = PetscViewerHDF5PopGroup(viewer); CHKERRQ(ierr);
    
    // then save explicit diffusion terms
    ierr = PetscViewerHDF5PushGroup(viewer, "diffusion"); CHKERRQ(ierr);
    for(unsigned int i=0; i<diffCoeffs->nExplicit; ++i)
    {
        ierr = VecView(diff[i], viewer); CHKERRQ(ierr);
    }
    ierr = PetscViewerHDF5PopGroup(viewer); CHKERRQ(ierr);
    
    // destroy viewer
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

    ierr = PetscLogStagePop(); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}


PetscErrorCode NavierStokesSolver::readRestartData(const std::string &filePath)
{
    PetscFunctionBeginUser;
    
    PetscErrorCode  ierr;
    
    PetscBool   fileExist = PETSC_FALSE;
    
    PetscViewer     viewer;

    // test if file exist 
    ierr = PetscTestFile((filePath+".h5").c_str(), 'w', &fileExist); CHKERRQ(ierr);
    
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
    
    // write extra data to the end of the file
    // save explicit convection terms first
    ierr = PetscViewerHDF5PushGroup(viewer, "convection"); CHKERRQ(ierr);
    for(unsigned int i=0; i<convCoeffs->nExplicit; ++i)
    {
        ierr = VecLoad(conv[i], viewer); CHKERRQ(ierr);
    }
    ierr = PetscViewerHDF5PopGroup(viewer); CHKERRQ(ierr);
    
    // then save explicit diffusion terms
    ierr = PetscViewerHDF5PushGroup(viewer, "diffusion"); CHKERRQ(ierr);
    for(unsigned int i=0; i<diffCoeffs->nExplicit; ++i)
    {
        ierr = VecLoad(diff[i], viewer); CHKERRQ(ierr);
    }
    ierr = PetscViewerHDF5PopGroup(viewer); CHKERRQ(ierr);
    
    // destroy viewer
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
    
    // update values and Eqs of ghost points based on current solution
    // TODO: for convective BCs, it's not totally correct
    ierr = bc->setGhostICs(solution); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}


// write numbers of iterations and residuals of linear solvers to a file
PetscErrorCode NavierStokesSolver::writeIterations(
        const int &timeIndex, const std::string &filePath)
{
    PetscErrorCode ierr;
    PetscInt nIters;
    PetscReal res;
    std::ofstream outfile;

    PetscFunctionBeginUser;

    if (mesh->mpiRank == 0)
    {
        if (timeIndex == 1)
            outfile.open(filePath.c_str());
        else
            outfile.open(filePath.c_str(), std::ios::out | std::ios::app);
        
        outfile << timeIndex << '\t';
        
        ierr = vSolver->getIters(nIters); CHKERRQ(ierr);
        ierr = vSolver->getResidual(res); CHKERRQ(ierr);
        outfile << nIters << '\t' << res << '\t';
        
        ierr = pSolver->getIters(nIters); CHKERRQ(ierr);
        ierr = pSolver->getResidual(res); CHKERRQ(ierr);
        outfile << nIters << '\t' << res << std::endl;
        
        outfile.close();
    }

    PetscFunctionReturn(0);
} // writeIterations


// manual finalization
PetscErrorCode NavierStokesSolver::finalize()
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
