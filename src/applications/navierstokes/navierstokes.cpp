/**
 * \file navierstokes.cpp
 * \brief Implementation of the class \c NavierStokesSolver.
 */

#include <fstream>
#include <string>

#include "navierstokes.h"


namespace petibm
{
namespace applications
{

NavierStokesSolver::NavierStokesSolver()
{

} // NavierStokesSolver


NavierStokesSolver::NavierStokesSolver(
		utilities::CartesianMesh &mesh2,
		utilities::FlowDescription &flow2,
		utilities::SimulationParameters &parameters2)
{
	mesh = mesh2;
	flow = flow2;
	parameters = parameters2;
	// register events
	PetscLogStageRegister("initialize", &stageInitialize);
	PetscLogStageRegister("rhsVelocity", &stageRHSVelocity);
	PetscLogStageRegister("solveVelocity", &stageSolveVelocity);
	PetscLogStageRegister("rhsPoisson", &stageRHSPoisson);
	PetscLogStageRegister("solvePoisson", &stageSolvePoisson);
	PetscLogStageRegister("projectionStep", &stageProjectionStep);
	PetscLogStageRegister("write", &stageWrite);
} // NavierStokesSolver


NavierStokesSolver::~NavierStokesSolver()
{

} // ~NavierStokesSolver


PetscErrorCode NavierStokesSolver::initialize()
{
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	ierr = PetscLogStagePush(stageInitialize); CHKERRQ(ierr);

	// set initial solution
	ierr = solution.init(mesh, parameters.output.format); CHKERRQ(ierr);

	// initialize boundary conditions
	ierr = bc.init(mesh); CHKERRQ(ierr);
	ierr = bc.setGhostICs(solution); CHKERRQ(ierr);

	// initialize time schemes
	ierr = convection.init(parameters.schemes.convection); CHKERRQ(ierr);
	ierr = diffusion.init(parameters.schemes.diffusion); CHKERRQ(ierr);

	ierr = assembleOperators(); CHKERRQ(ierr);

	// allocate PETSc Vecs
	ierr = VecDuplicate(solution.pGlobal, &phi); CHKERRQ(ierr);
	ierr = VecDuplicate(solution.UGlobal, &bc1); CHKERRQ(ierr);
	ierr = VecDuplicate(solution.UGlobal, &rhs1); CHKERRQ(ierr);
	ierr = VecDuplicate(solution.pGlobal, &rhs2); CHKERRQ(ierr);
	ierr = VecDuplicate(solution.UGlobal, &gradP); CHKERRQ(ierr);
	Conv.resize(convection.nExplicit);
	for (unsigned int i=0; i<Conv.size(); ++i) {
		ierr = VecDuplicate(solution.UGlobal, &Conv[i]); CHKERRQ(ierr);
	}
	Diff.resize(diffusion.nExplicit);
	for (unsigned int i=0; i<Diff.size(); ++i) {
		ierr = VecDuplicate(solution.UGlobal, &Diff[i]); CHKERRQ(ierr);
	}

	{
		MatNullSpace nsp;
		ierr = MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, nullptr,
		                          &nsp); CHKERRQ(ierr);
		ierr = MatSetNullSpace(DBNG, nsp); CHKERRQ(ierr);
		ierr = MatNullSpaceDestroy(&nsp); CHKERRQ(ierr);
	}

	ierr = linsolvers::createLinSolver("velocity", parameters.vSolver.config,
	                                   parameters.vSolver.type,
	                                   vSolver); CHKERRQ(ierr);
	ierr = vSolver->setMatrix(A); CHKERRQ(ierr);

	ierr = linsolvers::createLinSolver("poisson", parameters.pSolver.config,
	                                   parameters.pSolver.type,
	                                   pSolver); CHKERRQ(ierr);
	ierr = pSolver->setMatrix(DBNG); CHKERRQ(ierr);

	ierr = PetscLogStagePop();

	PetscFunctionReturn(0);
} // initialize


PetscErrorCode NavierStokesSolver::assembleOperators()
{
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	ierr = operators::createIdentity(mesh, I); CHKERRQ(ierr);
	ierr = operators::createR(mesh, R); CHKERRQ(ierr);
	ierr = operators::createRInv(mesh, RInv); CHKERRQ(ierr);
	ierr = operators::createM(mesh, M); CHKERRQ(ierr);
	ierr = operators::createMHead(mesh, MHat); CHKERRQ(ierr);
	ierr = operators::createLaplacian(mesh, bc, L, LCorrection); CHKERRQ(ierr);
	ierr = operators::createBnHead(L, parameters.step.dt,
	                               diffusion.implicitCoeff * flow.nu, 1,
	                               BNHat); CHKERRQ(ierr);
	ierr = operators::createDivergence(mesh, bc, D, DCorrection,
	                                   PETSC_FALSE); CHKERRQ(ierr);
	ierr = operators::createGradient(mesh, G, PETSC_FALSE); CHKERRQ(ierr);
	ierr = operators::createConvection(mesh, bc, N); CHKERRQ(ierr);

	ierr = MatMatMult(BNHat, G, MAT_INITIAL_MATRIX, PETSC_DEFAULT,
	                  &BNG); CHKERRQ(ierr);
	ierr = MatMatMult(D, BNG, MAT_INITIAL_MATRIX,
	                  PETSC_DEFAULT, &DBNG); CHKERRQ(ierr);
	ierr = MatScale(I, 1.0/parameters.step.dt); CHKERRQ(ierr);
	ierr = MatDuplicate(I, MAT_COPY_VALUES, &A); CHKERRQ(ierr);
	ierr = MatAXPY(A, - diffusion.implicitCoeff * flow.nu,
	               L, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // assembleOperators


PetscErrorCode NavierStokesSolver::solve()
{
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	ierr = assembleRHSVelocity(); CHKERRQ(ierr);
	ierr = solveVelocity(); CHKERRQ(ierr);
	ierr = assembleRHSPoisson(); CHKERRQ(ierr);
	ierr = solvePoisson(); CHKERRQ(ierr);
	ierr = projectionStep(); CHKERRQ(ierr);
	ierr = bc.updateGhostValues(solution); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // solve


PetscErrorCode NavierStokesSolver::assembleRHSVelocity()
{
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	ierr = PetscLogStagePush(stageRHSVelocity); CHKERRQ(ierr);

	// set explicit part of time derivative ($\frac{u^n}{\Delta t}$)
	ierr = VecCopy(solution.UGlobal, rhs1); CHKERRQ(ierr);
	ierr = VecScale(rhs1, 1.0 / parameters.step.dt); CHKERRQ(ierr);

	// prepare convection terms from time index n
	for(unsigned int i=Conv.size()-1; i>0; i--) {
		ierr = VecSwap(Conv[i], Conv[i-1]); CHKERRQ(ierr);
	}
	ierr = MatMult(N, solution.UGlobal, Conv[0]); CHKERRQ(ierr);
	ierr = VecScale(Conv[0], -1.0); CHKERRQ(ierr);
	// add all explicit convective terms to RHS
	for(unsigned int i=0; i<Conv.size(); ++i) {
		ierr = VecAXPY(rhs1, convection.explicitCoeffs[i], Conv[i]); CHKERRQ(ierr);
	}

	// prepare diffusion terms from time index n
	for(unsigned int i=Diff.size()-1; i>0; i--) {
		ierr = VecSwap(Diff[i], Diff[i-1]); CHKERRQ(ierr);
	}
	ierr = MatMult(L, solution.UGlobal, Diff[0]); CHKERRQ(ierr);
	ierr = MatMultAdd(LCorrection, solution.UGlobal, Diff[0],
	                  Diff[0]); CHKERRQ(ierr);
	ierr = VecScale(Diff[0], flow.nu); CHKERRQ(ierr);
	// add all explicit diffusive terms to RHS
	for(unsigned int i=0; i<Diff.size(); ++i) {
		ierr = VecAXPY(rhs1, diffusion.explicitCoeffs[i], Diff[i]); CHKERRQ(ierr);
	}

	// add explicit pressure gradient to RHS
	ierr = MatMult(G, solution.pGlobal, gradP); CHKERRQ(ierr);
	ierr = VecScale(gradP, -1.0); CHKERRQ(ierr);
	ierr = VecAXPY(rhs1, 1.0, gradP); CHKERRQ(ierr);

	// update the Eq.s of ghost points
	ierr = bc.updateEqs(solution, parameters.step.dt); CHKERRQ(ierr);
	// prepare the diffusion BC correction from time index n+1
	ierr = MatMult(LCorrection, solution.UGlobal, bc1); CHKERRQ(ierr);
	ierr = VecScale(bc1, diffusion.implicitCoeff * flow.nu); CHKERRQ(ierr);
	// add inhomogeneous boundary values to RHS
	ierr = VecAXPY(rhs1, 1.0, bc1); CHKERRQ(ierr);

	ierr = PetscLogStagePop(); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // assembleRHSVelocity


PetscErrorCode NavierStokesSolver::solveVelocity()
{
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	ierr = PetscLogStagePush(stageSolveVelocity); CHKERRQ(ierr);

	ierr = vSolver->solve(solution.UGlobal, rhs1); CHKERRQ(ierr);

	ierr = PetscLogStagePop(); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // solveVelocity


PetscErrorCode NavierStokesSolver::assembleRHSPoisson()
{
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	ierr = PetscLogStagePush(stageRHSPoisson); CHKERRQ(ierr);

	ierr = MatMult(D, solution.UGlobal, rhs2); CHKERRQ(ierr);
	ierr = MatMultAdd(DCorrection, solution.UGlobal, rhs2, rhs2); CHKERRQ(ierr);

	ierr = PetscLogStagePop(); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // assembleRHSPoisson


PetscErrorCode NavierStokesSolver::solvePoisson()
{
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	ierr = PetscLogStagePush(stageSolvePoisson); CHKERRQ(ierr);

	// solve for the pressure correction
	ierr = pSolver->solve(phi, rhs2); CHKERRQ(ierr);

	ierr = PetscLogStagePop(); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // solvePoisson


PetscErrorCode NavierStokesSolver::projectionStep()
{
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	ierr = PetscLogStagePush(stageProjectionStep); CHKERRQ(ierr);

	// project velocity field onto divergence-free space
	ierr = MatMult(BNG, phi, rhs1); CHKERRQ(ierr);
	ierr = VecAXPY(solution.UGlobal, -1.0, rhs1); CHKERRQ(ierr);
	// correct pressure field
	ierr = VecAXPY(solution.pGlobal, 1.0, phi); CHKERRQ(ierr);

	ierr = PetscLogStagePop(); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // projectionStep


PetscErrorCode NavierStokesSolver::write(std::string directory,
                                         std::string fileName)
{
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	ierr = PetscLogStagePush(stageWrite); CHKERRQ(ierr);

	ierr = solution.write(directory, fileName); CHKERRQ(ierr);

	ierr = PetscLogStagePop(); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // write


PetscErrorCode NavierStokesSolver::writeIterations(int timeIndex,
                                                   std::string filePath)
{
	PetscErrorCode ierr;
	PetscMPIInt rank;
	PetscInt nIters;
	std::ofstream outfile;

	PetscFunctionBeginUser;

	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

	if (rank == 0)
	{
		if (timeIndex == 1)
			outfile.open(filePath.c_str());
		else
			outfile.open(filePath.c_str(), std::ios::out | std::ios::app);
		outfile << timeIndex << '\t';
		ierr = vSolver->getIters(nIters); CHKERRQ(ierr);
		outfile << nIters << '\t';
		ierr = pSolver->getIters(nIters); CHKERRQ(ierr);
		outfile << nIters << std::endl;
		outfile.close();
	}

	PetscFunctionReturn(0);
} // writeIterations


PetscErrorCode NavierStokesSolver::finalize()
{
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	ierr = VecDestroy(&phi); CHKERRQ(ierr);
	ierr = VecDestroy(&bc1); CHKERRQ(ierr);
	ierr = VecDestroy(&rhs1); CHKERRQ(ierr);
	ierr = VecDestroy(&rhs2); CHKERRQ(ierr);
	ierr = VecDestroy(&gradP); CHKERRQ(ierr);
	for (unsigned int i=0; i<Conv.size(); ++i) {
		ierr = VecDestroy(&Conv[i]); CHKERRQ(ierr);
	}
	for (unsigned int i=0; i<Diff.size(); ++i) {
		ierr = VecDestroy(&Diff[i]); CHKERRQ(ierr);
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
	ierr = MatDestroy(&BNHat); CHKERRQ(ierr);
	ierr = MatDestroy(&MHat); CHKERRQ(ierr);
	ierr = MatDestroy(&M); CHKERRQ(ierr);
	ierr = MatDestroy(&RInv); CHKERRQ(ierr);
	ierr = MatDestroy(&R); CHKERRQ(ierr);
	ierr = MatDestroy(&I); CHKERRQ(ierr);

	vSolver.~shared_ptr();
	pSolver.~shared_ptr();

	PetscFunctionReturn(0);
} // finalize

} // end of namespace applications
} // end of namespace petibm
