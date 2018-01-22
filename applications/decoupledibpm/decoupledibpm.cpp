/**
 * \file decoupledibpm.cpp
 * \brief Implementation of the class \c DecoupledIBPMSolver.
 */

// PETSc
# include <petscviewerhdf5.h>

// Decoupled solver
# include "decoupledibpm.h"


DecoupledIBPMSolver::DecoupledIBPMSolver(
        const petibm::type::Mesh &inMesh,
        const petibm::type::Boundary &inBC,
        const petibm::type::BodyPack &inBodies,
        const YAML::Node &node)
{
    initialize(inMesh, inBC, inBodies, node);
} // DecoupledIBPMSolver


DecoupledIBPMSolver::~DecoupledIBPMSolver()
{
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    PetscBool finalized;

    ierr = PetscFinalized(&finalized); CHKERRV(ierr);
    if (finalized) return;

    ierr = VecDestroy(&df); CHKERRV(ierr);
    ierr = VecDestroy(&f); CHKERRV(ierr);
    ierr = VecDestroy(&Eu); CHKERRV(ierr);
    
    ierr = MatDestroy(&H); CHKERRV(ierr);
    ierr = MatDestroy(&E); CHKERRV(ierr);
    ierr = MatDestroy(&EBNH); CHKERRV(ierr);
    ierr = MatDestroy(&BNH); CHKERRV(ierr);
}


// destroy
PetscErrorCode DecoupledIBPMSolver::destroy()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    
    ierr = VecDestroy(&df); CHKERRQ(ierr);
    ierr = VecDestroy(&f); CHKERRQ(ierr);
    ierr = VecDestroy(&Eu); CHKERRQ(ierr);
    
    ierr = MatDestroy(&H); CHKERRQ(ierr);
    ierr = MatDestroy(&E); CHKERRQ(ierr);
    ierr = MatDestroy(&EBNH); CHKERRQ(ierr);
    ierr = MatDestroy(&BNH); CHKERRQ(ierr);
    
    fSolver.reset();
    bodies.reset();
    
    ierr = NavierStokesSolver::destroy(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // finalize


PetscErrorCode DecoupledIBPMSolver::initialize(
        const petibm::type::Mesh &inMesh,
        const petibm::type::Boundary &inBC,
        const petibm::type::BodyPack &inBodies,
        const YAML::Node &node) 
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    
    ierr = NavierStokesSolver::initialize(inMesh, inBC, node); CHKERRQ(ierr);
    
    // continue on stageInitialize
    ierr = PetscLogStagePush(stageInitialize); CHKERRQ(ierr);
    
    // make a reference to provided BodyPack
    bodies = inBodies;

    // create force solver
    ierr = petibm::linsolver::createLinSolver("forces", settings, fSolver); CHKERRQ(ierr);

    // create extra operators for decoupled method
    ierr = createExtraOperators(); CHKERRQ(ierr);

    // create extra vectors for decoupled method
    ierr = createExtraVectors(); CHKERRQ(ierr);

    // set coefficient matrix to the force solver
    ierr = fSolver->setMatrix(EBNH); CHKERRQ(ierr);
    
    // register additional events
    ierr = PetscLogStageRegister("rhsForces", &stageRHSForces); CHKERRQ(ierr);
    ierr = PetscLogStageRegister("solveForces", &stageSolveForces); CHKERRQ(ierr);
    ierr = PetscLogStageRegister("integrateForces", &stageIntegrateForces); CHKERRQ(ierr);
    
    ierr = PetscLogStagePop(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // initialize


// create extra operators for decoupled method
PetscErrorCode DecoupledIBPMSolver::createExtraOperators()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    
    Mat     R;
    Mat     MHat;
    Mat     BN;
    
    Vec     RDiag;
    Vec     MHatDiag;

    // create diagonal matrix R and get a Vec holding the diagonal
    ierr = petibm::operators::createR(mesh, R); CHKERRQ(ierr);
    ierr = MatCreateVecs(R, nullptr, &RDiag); CHKERRQ(ierr);
    ierr = MatGetDiagonal(R, RDiag); CHKERRQ(ierr);
    
    // create diagonal matrix M and get a Vec holding the diagonal
    ierr = petibm::operators::createMHead(mesh, MHat); CHKERRQ(ierr);
    ierr = MatCreateVecs(MHat, nullptr, &MHatDiag); CHKERRQ(ierr);
    ierr = MatGetDiagonal(MHat, MHatDiag); CHKERRQ(ierr);

    // create a Delta operator and its transpose
    ierr = petibm::operators::createDelta(mesh, bc, bodies, E); CHKERRQ(ierr);
    ierr = MatTranspose(E, MAT_INITIAL_MATRIX, &H); CHKERRQ(ierr);
    
    // create operator E
    ierr = MatDiagonalScale(E, nullptr, RDiag); CHKERRQ(ierr);
    ierr = MatDiagonalScale(E, nullptr, MHatDiag); CHKERRQ(ierr);
    
    // create operator H
    // Note: we do nothing here. Because if f is force (unit [m/s^2]), we need to 
    // convert it to pressure when using it in Eulerian momentum equations.
    // That is, (H*R^{-1}*f). While H = R*Delta, This means the real scattering
    // operator used in momentum equation is just a Delta operator. So we just
    // let H be the Delta operator. One thing to remember is that this is based
    // on the assumption that the size of Lagragian element is equal to the size
    // of the Eulerian grid near by.

    // create operator BN
    ierr = petibm::operators::createBnHead(
            L, dt, diffCoeffs->implicitCoeff*nu, 1, BN); CHKERRQ(ierr);

    // create BNH
    ierr = MatMatMult(
            BN, H, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &BNH); CHKERRQ(ierr);
    
    // create EBNH
    ierr = MatMatMult(
            E, BNH, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &EBNH); CHKERRQ(ierr);
    
    
    // destroy temporary Vecs and Mats
    ierr = VecDestroy(&RDiag); CHKERRQ(ierr);
    ierr = VecDestroy(&MHatDiag); CHKERRQ(ierr);
    ierr = MatDestroy(&MHat); CHKERRQ(ierr);
    ierr = MatDestroy(&R); CHKERRQ(ierr);
    ierr = MatDestroy(&BN); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // assembleOperators


// create extra vectors for decoupled method
PetscErrorCode DecoupledIBPMSolver::createExtraVectors()
{
    PetscFunctionBeginUser;
    
    PetscErrorCode  ierr;
    
    // create f and Eu 
    ierr = MatCreateVecs(EBNH, &f, &Eu); CHKERRQ(ierr);
    
    // create df
    ierr = VecDuplicate(f, &df); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}


// advance one time-step
PetscErrorCode DecoupledIBPMSolver::advance()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    // prepare velocity system and solve it
    ierr = assembleRHSVelocity(); CHKERRQ(ierr);
    ierr = solveVelocity(); CHKERRQ(ierr);
    
    // prepare force system and solve it
    ierr = assembleRHSForces(); CHKERRQ(ierr);
    ierr = solveForces(); CHKERRQ(ierr);
    
    // prepare poisson system and solve it
    ierr = assembleRHSPoisson(); CHKERRQ(ierr);
    ierr = solvePoisson(); CHKERRQ(ierr);
    
    // project solution to divergence free field
    ierr = projectionStep(); CHKERRQ(ierr);
    
    // update values of ghost points
    ierr = bc->updateGhostValues(solution); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // solve


// prepare right-hand-side for velocity system
PetscErrorCode DecoupledIBPMSolver::assembleRHSVelocity()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    // prepare the part coming from underlying Navier-Stoke solver
    ierr = NavierStokesSolver::assembleRHSVelocity(); CHKERRQ(ierr);

    // continue on stageRHSVelocity
    ierr = PetscLogStagePush(stageRHSVelocity); CHKERRQ(ierr);
    
    // prepare distributed boundary forces
    ierr = MatMultAdd(H, f, rhs1, rhs1); CHKERRQ(ierr);

    ierr = PetscLogStagePop(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // assembleRHSVelocity


// prepare right-hand-side for velocity system
PetscErrorCode DecoupledIBPMSolver::assembleRHSForces()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = PetscLogStagePush(stageRHSForces); CHKERRQ(ierr);

    // right-hand-side is -Eu^{**}
    ierr = MatMult(E, solution->UGlobal, Eu); CHKERRQ(ierr);
    ierr = VecScale(Eu, -1.0); CHKERRQ(ierr);
    
    // zerolize df, because the force increment should not have anything to do
    // with that of the previous teim-step
    ierr = VecSet(df, 0.0); CHKERRQ(ierr);

    ierr = PetscLogStagePop(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // assembleRHSForces


// solve force increment
PetscErrorCode DecoupledIBPMSolver::solveForces()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = PetscLogStagePush(stageSolveForces); CHKERRQ(ierr);

    // solve for the forces correction
    ierr = fSolver->solve(df, Eu); CHKERRQ(ierr);

    ierr = PetscLogStagePop(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // solveForces


// prepare right-hand-side for Poisson system
PetscErrorCode DecoupledIBPMSolver::assembleRHSPoisson()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    // continue on stageRHSPoisson
    ierr = PetscLogStagePush(stageRHSPoisson); CHKERRQ(ierr);
    
    // correct intermediate velocity with forcing term
    ierr = MatMultAdd(
            BNH, df, solution->UGlobal, solution->UGlobal); CHKERRQ(ierr);
    
    ierr = PetscLogStagePop(); CHKERRQ(ierr);

    // create rhs2 with the same way Navier-Stokes solver does
    ierr = NavierStokesSolver::assembleRHSPoisson(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // assembleRHSPoisson


// project solution to divergence-free field
PetscErrorCode DecoupledIBPMSolver::projectionStep()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;


    // project velocity field onto divergence-free space
    ierr = NavierStokesSolver::projectionStep(); CHKERRQ(ierr);
    
    // correct Lagrangian forces
    ierr = PetscLogStagePush(stageProjectionStep); CHKERRQ(ierr);
    
    ierr = VecAXPY(f, 1.0, df); CHKERRQ(ierr);
    
    ierr = PetscLogStagePop(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // projectionStep


// output extra data required for restarting to the user provided file
PetscErrorCode DecoupledIBPMSolver::writeRestartData(
  const PetscReal &t, const std::string &filePath)
{
    PetscFunctionBeginUser;
    
    PetscErrorCode      ierr;
    PetscViewer         viewer;
    
    ierr = NavierStokesSolver::writeRestartData(t, filePath); CHKERRQ(ierr);

    ierr = PetscLogStagePush(stageWrite); CHKERRQ(ierr);
    
    // create PetscViewer with append mode
    ierr = PetscViewerCreate(mesh->comm, &viewer); CHKERRQ(ierr);
    ierr = PetscViewerSetType(viewer, PETSCVIEWERHDF5); CHKERRQ(ierr);
    ierr = PetscViewerFileSetMode(viewer, FILE_MODE_APPEND); CHKERRQ(ierr);
    ierr = PetscViewerFileSetName(viewer, (filePath+".h5").c_str()); CHKERRQ(ierr);
    
    // go to the root node first. Not necessary. Just in case.
    ierr = PetscViewerHDF5PushGroup(viewer, "/"); CHKERRQ(ierr);
    
    // save Lagragian forces
    ierr = PetscObjectSetName((PetscObject) f, "force"); CHKERRQ(ierr);
    ierr = VecView(f, viewer); CHKERRQ(ierr);
    
    // destroy viewer
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

    ierr = PetscLogStagePop(); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}


// read data necessary for restarting
PetscErrorCode DecoupledIBPMSolver::readRestartData(
  const std::string &filePath, PetscReal &t)
{
    PetscFunctionBeginUser;
    
    PetscErrorCode      ierr;
    PetscViewer         viewer;
    
    ierr = NavierStokesSolver::readRestartData(filePath, t); CHKERRQ(ierr);
    
    // create PetscViewer with read-only mode
    ierr = PetscViewerCreate(mesh->comm, &viewer); CHKERRQ(ierr);
    ierr = PetscViewerSetType(viewer, PETSCVIEWERHDF5); CHKERRQ(ierr);
    ierr = PetscViewerFileSetMode(viewer, FILE_MODE_READ); CHKERRQ(ierr);
    ierr = PetscViewerFileSetName(viewer, (filePath+".h5").c_str()); CHKERRQ(ierr);
    
    // go to the root node first. Not necessary. Just in case.
    ierr = PetscViewerHDF5PushGroup(viewer, "/"); CHKERRQ(ierr);

    // save Lagragian forces
    ierr = PetscObjectSetName((PetscObject) f, "force"); CHKERRQ(ierr);
    ierr = VecLoad(f, viewer); CHKERRQ(ierr);
    
    // destroy viewer
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}


// write # of iterations and residuals of linear solvers
PetscErrorCode DecoupledIBPMSolver::writeIterations(
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
            asciiViewers[filePath], "\t%d\t%e", nIters, res); CHKERRQ(ierr);
    
    ierr = fSolver->getIters(nIters); CHKERRQ(ierr);
    ierr = fSolver->getResidual(res); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(
            asciiViewers[filePath], "\t%d\t%e\n", nIters, res); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
} // writeIterations


// write averaged forces
PetscErrorCode DecoupledIBPMSolver::writeIntegratedForces(
            const PetscReal &t, const std::string &filePath)
{
    PetscErrorCode ierr;
    petibm::type::RealVec2D     fAvg;

    PetscFunctionBeginUser;

    ierr = PetscLogStagePush(stageIntegrateForces); CHKERRQ(ierr);

    // get averaged forces first
    ierr = bodies->calculateAvgForces(f, fAvg); CHKERRQ(ierr);

    ierr = PetscLogStagePop(); CHKERRQ(ierr);
    
    // write current time
    ierr = PetscViewerASCIIPrintf(
            asciiViewers[filePath], "%10.8e", t); CHKERRQ(ierr);
    
    // write forces body by body
    for(int i=0; i<bodies->nBodies; ++i)
    {
        for(int d=0; d<mesh->dim; ++d)
        {
            ierr = PetscViewerASCIIPrintf(
                    asciiViewers[filePath], "\t%10.8e", fAvg[i][d]); CHKERRQ(ierr);
        }
    }
    ierr = PetscViewerASCIIPrintf(asciiViewers[filePath], "\n"); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
} // writeIntegratedForces
