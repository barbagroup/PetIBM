/*! Solve the flow around a 2D inline oscillating cylinder using an decoupled
 *  version of the immersed-boundary projection method (Li et al., 2016).
 * \file decoupledibpm.cpp
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 *
 * Li, R. Y., Xie, C. M., Huang, W. X., & Xu, C. X. (2016).
 * An efficient immersed boundary projection method for flow over
 * complex/moving boundaries.
 * Computers & Fluids, 140, 122-135.
 */

#include <iomanip>
#include <sys/stat.h>

#include <petscsys.h>
#include <yaml-cpp/yaml.h>

#include <petibm/parser.h>
#include <petibm/mesh.h>
#include <petibm/boundary.h>
#include <petibm/bodypack.h>
#include <petibm/solution.h>
#include <petibm/timeintegration.h>
#include <petibm/linsolver.h>
#include <petibm/operators.h>

/* \struct Kinematics
 * \brief A structure holding the parameters of the kinematics.
 */
struct Kinematics
{
	PetscReal Am;  ///< amplitude of the harmonic oscillation
	PetscReal f;   ///< oscillation frequency
	PetscReal D;   ///< characteristic length of the body
}; // Kinematics


/* \struct DecouplingParams
 * \brief A structure holding the parameters of the iterative procedure.
 */
struct DecouplingParams
{
	PetscInt maxIters = 1;     ///< Maximum number of inner iterations
	PetscReal atol = 1.0E-12;  ///< Criterion (absolute tolerance) to stop
}; // DecouplingParams


/* \fn AppParseBoundaryKinematics
 * \brief Parse the configuration file to get the parameters of the kinematics.
 *
 * \param filepath Path of the configuration file.
 * \param ctx Struct to hold the parameters (passed by reference).
 */
PetscErrorCode AppParseBoundaryKinematics(
	const std::string filepath, Kinematics &ctx)
{
	YAML::Node config;

	PetscFunctionBeginUser;

	config = YAML::LoadFile(filepath);
	if (config["bodies"].size() > 1)
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
		        "This example only works for a single immersed boundary");
	const YAML::Node &body = config["bodies"][0];
	ctx.f = body["frequency"].as<PetscReal>(0.2);
	ctx.Am = body["amplitude"].as<PetscReal>(1.0 / (2.0 * PETSC_PI * ctx.f));
	ctx.D = body["length"].as<PetscReal>(1.0);

	PetscFunctionReturn(0);
} // AppParseBoundaryKinematics


/* \fn AppParseDecouplingParams
 * \brief Parse the configuration file to get the parameters of the iterative
 *        procedure.
 *
 * \param filepath Path of the configuration file.
 * \param ctx Struct to hold the parameters (pass by reference).
 */
PetscErrorCode AppParseDecouplingParams(
	const std::string filepath, DecouplingParams &ctx)
{
	YAML::Node config;

	PetscFunctionBeginUser;

	config = YAML::LoadFile(filepath);
	const YAML::Node &params = config["parameters"]["decoupling"];
	ctx.maxIters = params["maxIters"].as<PetscInt>(1);
	ctx.atol = params["atol"].as<PetscReal>(1.0E-12);

	PetscFunctionReturn(0);
} // AppParseDecouplingParams


/* \fn AppSetBoundaryCoordinates
 * \brief Set the coordinates of the boundary at a given time.
 *
 * \param t The time.
 * \param ctx Struct with the parameters of the kinematics.
 * \param body The body which coordinates will be set (passed by reference).
 */
PetscErrorCode AppSetBoundaryCoordinates(
	const PetscReal t, const Kinematics ctx, petibm::type::SingleBody &body)
{
	PetscFunctionBeginUser;

	PetscReal Xc = -ctx.Am * PetscSinReal(2.0 * PETSC_PI * ctx.f * t),
	          Yc = 0.0;
	for (PetscInt k=0; k<body->nPts; k++)
	{
		body->coords[k][0] =
			Xc + 0.5 * ctx.D * PetscCosReal(2.0 * PETSC_PI * k / body->nPts);
		body->coords[k][1] =
			Yc + 0.5 * ctx.D * PetscSinReal(2.0 * PETSC_PI * k / body->nPts);
	}

	PetscFunctionReturn(0);
} // AppSetBoundaryCoordinates


/* \fn AppSetBoundaryVelocity
 * \brief Set the velocity of the boundary at a given time.
 *
 * \param t The time.
 * \param ctx Struct with the parameters of the kinematics.
 * \param body The immersed body.
 * \param Ub Vector for the boundary velocity.
 */
PetscErrorCode AppSetBoundaryVelocity(
	const PetscReal t, const Kinematics ctx, const petibm::type::SingleBody body,
	Vec Ub)
{
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	PetscReal Um = ctx.Am * 2.0 * PETSC_PI * ctx.f;  // maximum velocity
	PetscReal Ubx = -Um * PetscCosReal(2.0 * PETSC_PI * ctx.f * t);
	PetscReal **Ub_arr;
	ierr = DMDAVecGetArrayDOF(body->da, Ub, &Ub_arr); CHKERRQ(ierr);
	for (PetscInt k=body->bgPt; k<body->edPt; k++)
		Ub_arr[k][0] = Ubx;
	ierr = DMDAVecRestoreArrayDOF(body->da, Ub, &Ub_arr); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // AppSetBoundaryVelocity


/* \fn AppWriteBoundaryCoordinates
 * \brief Write the boundary coordinates into a file in ASCII format.
 *
 * \param filepath Path of the output file.
 * \param body The immersed body.
 */
PetscErrorCode AppWriteBoundaryCoordinates(
	const std::string filepath, const petibm::type::SingleBody body)
{
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	PetscViewer viewer;
	ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr);
	ierr = PetscViewerSetType(viewer, PETSCVIEWERASCII); CHKERRQ(ierr);
	ierr = PetscViewerFileSetMode(viewer, FILE_MODE_WRITE); CHKERRQ(ierr);
	ierr = PetscViewerFileSetName(viewer, filepath.c_str()); CHKERRQ(ierr);
	for (PetscInt k=0; k<body->nPts; k++)
	{
		ierr = PetscViewerASCIIPrintf(viewer, "%10.8e\t%10.8e\n",
		                              body->coords[k][0],
		                              body->coords[k][1]); CHKERRQ(ierr);
	}
	ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

	PetscFunctionReturn(0);
} // AppWriteBoundaryCoordinates


/*! Main function: compute the flow around an inline-oscillating cylinder using
 *  a decoupled immersed boundary projection method (Li et al., 2016).
 */
int main(int argc, char **argv)
{
	PetscErrorCode ierr;
	Kinematics kinematics;
	DecouplingParams decoupling;
	YAML::Node config;
	std::string directory;
	petibm::type::Mesh mesh;
	petibm::type::Boundary boundary;
	petibm::type::Solution solution;
	petibm::type::TimeIntegration convCoeffs, diffCoeffs;
	Mat D, DCorrection, G, L, LCorrection, A, BN, BNG, DBNG, N;
	petibm::type::LinSolver vSolver, pSolver, fSolver;
	PetscBool isRefP = PETSC_FALSE;
	std::vector<Vec> conv, diff;
	Vec dp;
	Vec rhs1, bc1, rhs2;
	PetscViewer viewerIters, viewerForces, viewerLog;
	petibm::type::BodyPack bodyPack;
	Vec f, df, rhsf, UB;
	Mat E = PETSC_NULL, H = PETSC_NULL, BNH = PETSC_NULL, EBNH = PETSC_NULL;

	PetscLogStage stageInitialize,
	              stageRHSVelocity,
	              stageSolveVelocity,
	              stageRHSForces,
	              stageSolveForces,
	              stageIntegrateForces,
	              stageRHSPoisson,
	              stageSolvePoisson,
	              stageUpdate,
	              stageMoveBoundary,
	              stageWrite;
	
	ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);

	ierr = PetscPrintf(
		PETSC_COMM_WORLD,
		"\n\n*** Inline-oscillating body using the decoupled IBPM ***\n\n");
	CHKERRQ(ierr);

	ierr = PetscLogStageRegister(
		"initialize", &stageInitialize); CHKERRQ(ierr);
	ierr = PetscLogStagePush(stageInitialize); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD, "Initializing ...\n"); CHKERRQ(ierr);

	// Parse input YAML configuration file
	ierr = petibm::parser::getSettings(config); CHKERRQ(ierr);
	directory = config["directory"].as<std::string>();
	const YAML::Node &params = config["parameters"];
	PetscInt nstart = params["startStep"].as<PetscInt>(0);
	PetscInt nt = params["nt"].as<PetscInt>();
	PetscInt nsave = params["nsave"].as<PetscInt>();
	PetscReal dt = params["dt"].as<PetscReal>();
	PetscReal t = nstart * dt;
	PetscReal nu = config["flow"]["nu"].as<PetscReal>();

	// Parse configuration file to get parameters of the kinematics
	// and of the iterative procedure
	ierr = AppParseDecouplingParams(
		config["config.yaml"].as<std::string>(), decoupling); CHKERRQ(ierr);
	ierr = AppParseBoundaryKinematics(
		config["config.yaml"].as<std::string>(), kinematics); CHKERRQ(ierr);

	// Create a Cartesian mesh
	ierr = petibm::mesh::createMesh(
		PETSC_COMM_WORLD, config, mesh); CHKERRQ(ierr);
	ierr = mesh->write(directory + "/grid"); CHKERRQ(ierr);
	
	// Set up the external boundary conditions and apply initial conditions
	ierr = petibm::boundary::createBoundary(
		mesh, config, boundary); CHKERRQ(ierr);
	ierr = petibm::solution::createSolution(mesh, solution); CHKERRQ(ierr);
	ierr = solution->applyIC(config); CHKERRQ(ierr);
	ierr = boundary->setGhostICs(solution); CHKERRQ(ierr);

	// Get the coefficients of the time-integration scheme
	ierr = petibm::timeintegration::createTimeIntegration(
		"convection", config, convCoeffs); CHKERRQ(ierr);
	ierr = petibm::timeintegration::createTimeIntegration(
		"diffusion", config, diffCoeffs); CHKERRQ(ierr);
	
	// Assemble the divergence operator and its boundary correction
	ierr = petibm::operators::createDivergence(
		mesh, boundary, D, DCorrection, PETSC_FALSE); CHKERRQ(ierr);
	ierr = PetscObjectRegisterDestroy((PetscObject) D); CHKERRQ(ierr);
	ierr = PetscObjectRegisterDestroy((PetscObject) DCorrection); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions(
		(PetscObject) D, nullptr, "-D_mat_view"); CHKERRQ(ierr);

	// Assemble the gradient operator
	ierr = petibm::operators::createGradient(
		mesh, G, PETSC_FALSE); CHKERRQ(ierr);
	ierr = PetscObjectRegisterDestroy((PetscObject) G); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions(
		(PetscObject) G, nullptr, "-G_mat_view"); CHKERRQ(ierr);

	// Assemble the Laplacian operator and its boundary correction
	ierr = petibm::operators::createLaplacian(
		mesh, boundary, L, LCorrection); CHKERRQ(ierr);
	ierr = PetscObjectRegisterDestroy((PetscObject) L); CHKERRQ(ierr);
	ierr = PetscObjectRegisterDestroy((PetscObject) LCorrection); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions(
		(PetscObject) L, nullptr, "-L_mat_view"); CHKERRQ(ierr);

	// Assemble the implicit operator
	ierr = MatDuplicate(L, MAT_COPY_VALUES, &A); CHKERRQ(ierr);
	ierr = MatScale(A, -diffCoeffs->implicitCoeff * nu); CHKERRQ(ierr);
	ierr = MatShift(A, 1.0/dt); CHKERRQ(ierr);
	ierr = PetscObjectRegisterDestroy((PetscObject) A); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions(
		(PetscObject) A, nullptr, "-A_mat_view"); CHKERRQ(ierr);

	// Create linear solver for the velocity system
	ierr = petibm::linsolver::createLinSolver(
		"velocity", config, vSolver); CHKERRQ(ierr);
	ierr = vSolver->setMatrix(A); CHKERRQ(ierr);

	// Assemble the first-order approximation of the inverse
	// of the implicit operator
	ierr = petibm::operators::createIdentity(mesh, BN); CHKERRQ(ierr);
	ierr = MatScale(BN, dt); CHKERRQ(ierr);
	ierr = PetscObjectRegisterDestroy((PetscObject) BN); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions(
		(PetscObject) BN, nullptr, "-BN_mat_view"); CHKERRQ(ierr);

	// Assemble the operator to project the velocity of the divergence-free space
	ierr = MatMatMult(
		BN, G, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &BNG); CHKERRQ(ierr);
	ierr = PetscObjectRegisterDestroy((PetscObject) BNG); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions(
		(PetscObject) BNG, nullptr, "-BNG_mat_view"); CHKERRQ(ierr);

	// Assemble the Poisson operator
	ierr = MatMatMult(
		D, BNG, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &DBNG); CHKERRQ(ierr);
	ierr = PetscObjectRegisterDestroy((PetscObject) DBNG); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions(
		(PetscObject) DBNG, nullptr, "-DBNG_mat_view"); CHKERRQ(ierr);

	// Create a linear solver for the pressure
	ierr = petibm::linsolver::createLinSolver(
		"poisson", config, pSolver); CHKERRQ(ierr);
	// Set the nullspace for Poisson system with Neumann conditions
	{
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
			// Pin pressure at bottom-left corner of the domain
			PetscInt row[1] = {0};
			ierr = MatZeroRowsColumns(
				DBNG, 1, row, 1.0, nullptr, nullptr); CHKERRQ(ierr);
			isRefP = PETSC_TRUE;
		}
	}
	// Set the operator of the linear solver for the pressure
  ierr = pSolver->setMatrix(DBNG); CHKERRQ(ierr);

  // Assemble the non-linear operator
	ierr = petibm::operators::createConvection(mesh, boundary, N); CHKERRQ(ierr);
	ierr = PetscObjectRegisterDestroy((PetscObject) N); CHKERRQ(ierr);
	ierr = PetscObjectViewFromOptions(
		(PetscObject) N, nullptr, "-N_mat_view"); CHKERRQ(ierr);

	// Create the vectors to hold the convective and diffusive terms
	conv.resize(convCoeffs->nExplicit);
	for (unsigned int i=0; i<conv.size(); ++i)
	{
		ierr = VecDuplicate(solution->UGlobal, &conv[i]); CHKERRQ(ierr);
		ierr = PetscObjectRegisterDestroy((PetscObject) conv[i]); CHKERRQ(ierr);
	}
	diff.resize(diffCoeffs->nExplicit);
	for (unsigned int i=0; i<diff.size(); ++i)
	{
		ierr = VecDuplicate(solution->UGlobal, &diff[i]); CHKERRQ(ierr);
		ierr = PetscObjectRegisterDestroy((PetscObject) diff[i]); CHKERRQ(ierr);
	}
	// Create other vectors
	ierr = VecDuplicate(solution->pGlobal, &dp); CHKERRQ(ierr);
	ierr = PetscObjectRegisterDestroy((PetscObject) dp); CHKERRQ(ierr);
	ierr = VecDuplicate(solution->UGlobal, &rhs1); CHKERRQ(ierr);
	ierr = PetscObjectRegisterDestroy((PetscObject) rhs1); CHKERRQ(ierr);
	ierr = VecDuplicate(solution->UGlobal, &bc1); CHKERRQ(ierr);
	ierr = PetscObjectRegisterDestroy((PetscObject) bc1); CHKERRQ(ierr);
	ierr = VecDuplicate(solution->pGlobal, &rhs2); CHKERRQ(ierr);
	ierr = PetscObjectRegisterDestroy((PetscObject) rhs2); CHKERRQ(ierr);
	
	// Create a viewer to output the number of iterations for each solvers
	ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewerIters); CHKERRQ(ierr);
	ierr = PetscViewerSetType(viewerIters, PETSCVIEWERASCII); CHKERRQ(ierr);
	ierr = PetscViewerFileSetMode(viewerIters, FILE_MODE_WRITE); CHKERRQ(ierr);
	ierr = PetscViewerFileSetName(
		viewerIters, (directory+"/iterations.txt").c_str()); CHKERRQ(ierr);
	ierr = PetscObjectRegisterDestroy((PetscObject) viewerIters); CHKERRQ(ierr);
	// Create a viewer to output the forces acting on the immersed body
	ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewerForces); CHKERRQ(ierr);
	ierr = PetscViewerSetType(viewerForces, PETSCVIEWERASCII); CHKERRQ(ierr);
	ierr = PetscViewerFileSetMode(viewerForces, FILE_MODE_WRITE); CHKERRQ(ierr);
	ierr = PetscViewerFileSetName(
		viewerForces, (directory+"/forces.txt").c_str()); CHKERRQ(ierr);
	ierr = PetscObjectRegisterDestroy((PetscObject) viewerForces); CHKERRQ(ierr);
	// Create a viewer to save a summery of the PETSc logging
	// A logging file will be saved every time the field solution is saved
	ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewerLog); CHKERRQ(ierr);
	ierr = PetscViewerSetType(viewerLog, PETSCVIEWERASCII); CHKERRQ(ierr);
	ierr = PetscViewerFileSetMode(viewerLog, FILE_MODE_WRITE); CHKERRQ(ierr);
	ierr = PetscObjectRegisterDestroy((PetscObject) viewerLog); CHKERRQ(ierr);

	// Create the linear solver for the forcing system
	ierr = petibm::linsolver::createLinSolver(
		"forces", config, fSolver); CHKERRQ(ierr);
	// Create the immersed body and required vectors
	ierr = petibm::body::createBodyPack(mesh, config, bodyPack); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(bodyPack->dmPack, &f); CHKERRQ(ierr);
	ierr = PetscObjectRegisterDestroy((PetscObject) f); CHKERRQ(ierr);
	ierr = VecDuplicate(f, &df); CHKERRQ(ierr);
	ierr = PetscObjectRegisterDestroy((PetscObject) df); CHKERRQ(ierr);
	ierr = VecDuplicate(f, &UB); CHKERRQ(ierr);
	ierr = PetscObjectRegisterDestroy((PetscObject) UB); CHKERRQ(ierr);
	ierr = VecDuplicate(f, &rhsf); CHKERRQ(ierr);
	ierr = PetscObjectRegisterDestroy((PetscObject) rhsf); CHKERRQ(ierr);
	ierr = VecDuplicate(f, &UB); CHKERRQ(ierr);
	ierr = PetscObjectRegisterDestroy((PetscObject) UB); CHKERRQ(ierr);
	petibm::type::SingleBody &body = (bodyPack->bodies)[0];
	// Set the boundary coordinates and the boundary velocity
	ierr = AppSetBoundaryCoordinates(t, kinematics, body); CHKERRQ(ierr);
	ierr = AppSetBoundaryVelocity(t, kinematics, body, UB); CHKERRQ(ierr);

	ierr = PetscPrintf(
		PETSC_COMM_WORLD, "Initialization complete!\n"); CHKERRQ(ierr);
	ierr = PetscLogStagePop(); CHKERRQ(ierr); // End of stageInitialize

	// Write to file the initial field solution and boundary coordinates
	ierr = PetscLogStageRegister(
		"write", &stageWrite); CHKERRQ(ierr);
	ierr = PetscLogStagePush(stageWrite); CHKERRQ(ierr);
	ierr = PetscPrintf(
		PETSC_COMM_WORLD, "[time-step 0] Writing solution ... "); CHKERRQ(ierr);
	ierr = solution->write(directory + "/solution/0000000"); CHKERRQ(ierr);
	ierr = AppWriteBoundaryCoordinates(
		directory + "/solution/0000000.curve", body); CHKERRQ(ierr);
	ierr = PetscPrintf(
		PETSC_COMM_WORLD, "done\n"); CHKERRQ(ierr);
	ierr = PetscLogStagePop(); CHKERRQ(ierr); // End of stageWrite

	ierr = PetscLogStageRegister(
		"set velocity RHS", &stageRHSVelocity); CHKERRQ(ierr);
	ierr = PetscLogStageRegister(
		"solve velocity", &stageSolveVelocity); CHKERRQ(ierr);
	ierr = PetscLogStageRegister(
		"set forces RHS", &stageRHSForces); CHKERRQ(ierr);
	ierr = PetscLogStageRegister(
		"solve forces", &stageSolveForces); CHKERRQ(ierr);
	ierr = PetscLogStageRegister(
		"set Poisson RHS", &stageRHSPoisson); CHKERRQ(ierr);
	ierr = PetscLogStageRegister(
		"solve Poisson", &stageSolvePoisson); CHKERRQ(ierr);
	ierr = PetscLogStageRegister(
		"update", &stageUpdate); CHKERRQ(ierr);
	ierr = PetscLogStageRegister(
		"integrate forces", &stageIntegrateForces); CHKERRQ(ierr);
	ierr = PetscLogStageRegister(
		"move boundary", & stageMoveBoundary); CHKERRQ(ierr);

	// Start time integration
	for (PetscInt ite=nstart+1; ite<=nstart+nt; ite++)
	{
		t += dt;

		// Move the immersed boundary (update coordinates and velocity)
		ierr = PetscLogStagePush(stageMoveBoundary); CHKERRQ(ierr);
		ierr = AppSetBoundaryCoordinates(t, kinematics, body); CHKERRQ(ierr);
		ierr = AppSetBoundaryVelocity(t, kinematics, body, UB); CHKERRQ(ierr);
		ierr = body->updateMeshIdx(); CHKERRQ(ierr);
		// Re-compute operators
		if (E != PETSC_NULL) {ierr = MatDestroy(&E); CHKERRQ(ierr);}
		if (H != PETSC_NULL) {ierr = MatDestroy(&H); CHKERRQ(ierr);}
		if (BNH != PETSC_NULL) {ierr = MatDestroy(&BNH); CHKERRQ(ierr);}
		if (EBNH != PETSC_NULL) {ierr = MatDestroy(&EBNH); CHKERRQ(ierr);}
		ierr = petibm::operators::createDelta(
			mesh, boundary, bodyPack, E); CHKERRQ(ierr);
		ierr = MatTranspose(E, MAT_INITIAL_MATRIX, &H); CHKERRQ(ierr);
		{
			Mat R, MHat;
			ierr = petibm::operators::createR(mesh, R); CHKERRQ(ierr);
			ierr = petibm::operators::createMHead(mesh, MHat); CHKERRQ(ierr);
			Vec dR, dMHat;
			ierr = MatCreateVecs(R, nullptr, &dR); CHKERRQ(ierr);
			ierr = MatGetDiagonal(R, dR); CHKERRQ(ierr);
			ierr = MatDestroy(&R); CHKERRQ(ierr);
			ierr = MatCreateVecs(MHat, nullptr, &dMHat); CHKERRQ(ierr);
			ierr = MatGetDiagonal(MHat, dMHat); CHKERRQ(ierr);
			ierr = MatDestroy(&MHat); CHKERRQ(ierr);
			ierr = MatDiagonalScale(E, nullptr, dR); CHKERRQ(ierr);
			ierr = MatDiagonalScale(E, nullptr, dMHat); CHKERRQ(ierr);
			ierr = VecDestroy(&dR); CHKERRQ(ierr);
			ierr = VecDestroy(&dMHat); CHKERRQ(ierr);
		}
		ierr = MatMatMult(
			BN, H, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &BNH); CHKERRQ(ierr);
		ierr = MatMatMult(
			E, BNH, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &EBNH); CHKERRQ(ierr);
		ierr = PetscObjectViewFromOptions(
			(PetscObject) EBNH, nullptr, "-EBNH_mat_view"); CHKERRQ(ierr);
		ierr = fSolver->setMatrix(EBNH); CHKERRQ(ierr);
		ierr = PetscLogStagePop(); CHKERRQ(ierr); // End of stageMoveBoundary

		// Solve the forcing system to obtain prediction of the Lagrangian forces
		ierr = PetscLogStagePush(stageRHSForces); CHKERRQ(ierr);
		ierr = MatMult(E, solution->UGlobal, rhsf); CHKERRQ(ierr);
		ierr = VecAYPX(rhsf, -1.0, UB); CHKERRQ(ierr);
		ierr = PetscLogStagePop(); CHKERRQ(ierr); // End of stageRHSForces
		ierr = PetscLogStagePush(stageSolveForces); CHKERRQ(ierr);
		ierr = fSolver->solve(f, rhsf); CHKERRQ(ierr);
		ierr = PetscLogStagePop(); CHKERRQ(ierr); // End of stageSolveForces

		// Set up the right-hand side of the velocity system
		ierr = PetscLogStagePush(stageRHSVelocity); CHKERRQ(ierr);
		ierr = MatMult(G, solution->pGlobal, rhs1); CHKERRQ(ierr);
		ierr = VecScale(rhs1, -1.0); CHKERRQ(ierr);
		ierr = VecAXPY(rhs1, 1.0/dt, solution->UGlobal); CHKERRQ(ierr);
		for (unsigned int i=conv.size()-1; i>0; i--)
		{
			ierr = VecSwap(conv[i], conv[i-1]); CHKERRQ(ierr);
		}
		ierr = MatMult(N, solution->UGlobal, conv[0]); CHKERRQ(ierr);
		for (unsigned int i=0; i<conv.size(); i++)
		{
			ierr = VecAXPY(
				rhs1, -convCoeffs->explicitCoeffs[i], conv[i]); CHKERRQ(ierr);
		}
		ierr = MatMult(L, solution->UGlobal, diff[0]); CHKERRQ(ierr);
		ierr = MatMultAdd(
			LCorrection, solution->UGlobal, diff[0], diff[0]); CHKERRQ(ierr);
		ierr = VecScale(diff[0], nu); CHKERRQ(ierr);
		for (unsigned int i=0; i<diff.size(); i++)
		{
			ierr = VecAXPY(
				rhs1, diffCoeffs->explicitCoeffs[i], diff[i]); CHKERRQ(ierr);
		}
		ierr = boundary->updateEqs(solution, dt); CHKERRQ(ierr);
		ierr = MatMult(LCorrection, solution->UGlobal, bc1); CHKERRQ(ierr);
		ierr = VecScale(bc1, nu); CHKERRQ(ierr);
		ierr = VecAXPY(rhs1, diffCoeffs->implicitCoeff, bc1); CHKERRQ(ierr);
		ierr = MatMultAdd(H, f, rhs1, rhs1); CHKERRQ(ierr);
		ierr = PetscLogStagePop(); // End of stageRHSVelocity

		// Reset the forces and pressure corrections
		ierr = VecSet(df, 0.0); CHKERRQ(ierr);
		ierr = VecSet(dp, 0.0); CHKERRQ(ierr);

		// Start the inner iterative procedure
		PetscInt inner = 0;
		PetscReal norm = decoupling.atol + 1.0;
		PetscInt nIters_v = 0, nIters_f = 0;
		PetscReal res_v = 0.0, res_f = 0.0;
		while (norm > decoupling.atol and inner < decoupling.maxIters)
		{
			inner++;
			
			// Add explicit forces correction to the RHS of the velocity system
			if (inner > 1)
			{
				ierr = PetscLogStagePush(stageRHSVelocity); CHKERRQ(ierr);
				ierr = MatMultAdd(H, df, rhs1, rhs1); CHKERRQ(ierr);
				ierr = PetscLogStagePop(); CHKERRQ(ierr); // End of stageRHSVelocity
			}
	
			// Solve the system for the velocity
			ierr = PetscLogStagePush(stageSolveVelocity); CHKERRQ(ierr);
			ierr = vSolver->solve(solution->UGlobal, rhs1); CHKERRQ(ierr);
			ierr = PetscLogStagePop(); CHKERRQ(ierr); // End of stageSolveVelocity
	
			// Set up the RHS for the forcing system
			ierr = PetscLogStagePush(stageRHSForces); CHKERRQ(ierr);
			ierr = MatMult(E, solution->UGlobal, rhsf); CHKERRQ(ierr);
			ierr = VecAYPX(rhsf, -1.0, UB); CHKERRQ(ierr);
			ierr = PetscLogStagePop(); CHKERRQ(ierr); // End of stageRHSForces
	
			// Solve the system for the forces correction
			ierr = PetscLogStagePush(stageSolveForces); CHKERRQ(ierr);
			ierr = fSolver->solve(df, rhsf); CHKERRQ(ierr);
			ierr = PetscLogStagePop(); CHKERRQ(ierr); // End of stageSolveForces
	
			// Apply the no-slip condition onto the velocity field
			// and update the Lagrangian forces
			ierr = PetscLogStagePush(stageUpdate); CHKERRQ(ierr);
			ierr = MatMultAdd(
				BNH, df, solution->UGlobal, solution->UGlobal); CHKERRQ(ierr);
			ierr = VecAXPY(f, 1.0, df); CHKERRQ(ierr);
			ierr = PetscLogStagePop(); CHKERRQ(ierr); // End of stageUpdate

			// Hold the number of iterations and the residuals for output later
			{
				PetscInt nIters_tmp;
				PetscReal res_tmp;
				ierr = vSolver->getIters(nIters_tmp); CHKERRQ(ierr);
				ierr = vSolver->getResidual(res_tmp); CHKERRQ(ierr);
				nIters_v += nIters_tmp;
				res_v += res_tmp;
				ierr = fSolver->getIters(nIters_tmp); CHKERRQ(ierr);
				ierr = fSolver->getResidual(res_tmp); CHKERRQ(ierr);
				nIters_f += nIters_tmp;
				res_f += res_tmp;
			}

			// Check the L-inf norm of the forces correction
			ierr = VecNorm(df, NORM_INFINITY, &norm); CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_WORLD,
			                   "[time-step %d][inner %d] L_inf(df) = %10.8e\n",
			                   ite, inner, norm); CHKERRQ(ierr);
		}
		// Average number of iterations and residuals over number of sub-steps
		nIters_v /= inner; res_v /= inner;
		nIters_f /= inner; res_f /= inner;

		// Set up the RHS of the Poisson system
		ierr = PetscLogStagePush(stageRHSPoisson); CHKERRQ(ierr);
		ierr = MatMult(D, solution->UGlobal, rhs2); CHKERRQ(ierr);
		ierr = MatMultAdd(
			DCorrection, solution->UGlobal, rhs2, rhs2); CHKERRQ(ierr);
		if (isRefP) // Case where we pin the pressure at bottom-left corner
		{
			ierr = VecSetValue(rhs2, 0, 0.0, INSERT_VALUES); CHKERRQ(ierr);
			ierr = VecAssemblyBegin(rhs2); CHKERRQ(ierr);
			ierr = VecAssemblyEnd(rhs2); CHKERRQ(ierr);
		}
		ierr = PetscLogStagePop(); CHKERRQ(ierr); // End of stageRHSPoisson

		// Solver the Poisson system for the pressure
		ierr = PetscLogStagePush(stageSolvePoisson); CHKERRQ(ierr);
		ierr = pSolver->solve(dp, rhs2); CHKERRQ(ierr);
		ierr = PetscLogStagePop(); CHKERRQ(ierr); // End of stageSolvePoisson

		// Get the number of iterations and residuals for the Poisson solver
		PetscInt nIters_p;
		PetscReal res_p;
		ierr = pSolver->getIters(nIters_p); CHKERRQ(ierr);
		ierr = pSolver->getResidual(res_p); CHKERRQ(ierr);

		// Project the velocity field onto the divergence-free space
		// and update the pressure field
		ierr = PetscLogStagePush(stageUpdate); CHKERRQ(ierr);
		ierr = VecScale(solution->UGlobal, -1.0); CHKERRQ(ierr);
		ierr = MatMultAdd(
			BNG, dp, solution->UGlobal, solution->UGlobal); CHKERRQ(ierr);
		ierr = VecScale(solution->UGlobal, -1.0); CHKERRQ(ierr);
		ierr = VecAXPY(solution->pGlobal, 1.0, dp); CHKERRQ(ierr);
		ierr = PetscLogStagePop(); // End of stageUpdate

		// Write number of iterations and residuals into file
		ierr = PetscLogStagePush(stageWrite); CHKERRQ(ierr);
		ierr = PetscViewerASCIIPrintf(
			viewerIters, "%d\t%d\t", ite, inner); CHKERRQ(ierr);
		ierr = PetscViewerASCIIPrintf(
			viewerIters, "%d\t%e\t", nIters_v / inner, res_v / inner); CHKERRQ(ierr);
		ierr = PetscViewerASCIIPrintf(
			viewerIters, "%d\t%e\t", nIters_p, res_p); CHKERRQ(ierr);
		ierr = PetscViewerASCIIPrintf(
			viewerIters, "%d\t%e\n", nIters_f, res_f); CHKERRQ(ierr);
		ierr = PetscViewerFileSetMode(
			viewerIters, FILE_MODE_APPEND); CHKERRQ(ierr);
		ierr = PetscLogStagePop(); CHKERRQ(ierr); // End of stageWrite

		if (ite % nsave == 0)
		{
			// Write the Eulerian field solution and the body coordinates into files
			ierr = PetscLogStagePush(stageWrite); CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_WORLD,
			                   "[time-step %d] Writing solution ... ",
			                   ite); CHKERRQ(ierr);
			std::stringstream ss;
			ss << "/solution/" << std::setfill('0') << std::setw(7) << ite;
			ierr = solution->write(directory + ss.str()); CHKERRQ(ierr);
			ierr = AppWriteBoundaryCoordinates(
				directory + ss.str() + ".curve", body); CHKERRQ(ierr);
			ierr = PetscLogStagePop(); CHKERRQ(ierr); // End of stageWrite
			ierr = PetscViewerFileSetName(
				viewerLog, (directory + ss.str() + ".log").c_str()); CHKERRQ(ierr);
			ierr = PetscLogView(viewerLog); CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_WORLD, "done\n"); CHKERRQ(ierr);
		}

		// Compute the directional forces acting of the immersed body
		// and write them into file
		ierr = PetscLogStagePush(stageIntegrateForces); CHKERRQ(ierr);
		petibm::type::RealVec2D F;
		ierr = bodyPack->calculateAvgForces(f, F); CHKERRQ(ierr);
		ierr = PetscViewerASCIIPrintf(viewerForces, "%10.8e", t); CHKERRQ(ierr);
		for(int d=0; d<mesh->dim; d++)
		{
			ierr = PetscViewerASCIIPrintf(
				viewerForces, "\t%10.8e", F[0][d]); CHKERRQ(ierr);
		}
		ierr = PetscViewerASCIIPrintf(viewerForces, "\n"); CHKERRQ(ierr);
		ierr = PetscViewerFileSetMode(
			viewerForces, FILE_MODE_APPEND); CHKERRQ(ierr);
		ierr = PetscLogStagePop(); CHKERRQ(ierr); // End of stageIntegrateForces
	}

	// Free space
	ierr = MatDestroy(&E); CHKERRQ(ierr);
	ierr = MatDestroy(&H); CHKERRQ(ierr);
	ierr = MatDestroy(&BNH); CHKERRQ(ierr);
	ierr = MatDestroy(&EBNH); CHKERRQ(ierr);
	config.~Node();
	mesh.~shared_ptr();
	boundary.~shared_ptr();
	bodyPack.~shared_ptr();
	solution.~shared_ptr();
	convCoeffs.~shared_ptr();
	diffCoeffs.~shared_ptr();
	vSolver.~shared_ptr();
	pSolver.~shared_ptr();
	fSolver.~shared_ptr();

	ierr = PetscPrintf(
		PETSC_COMM_WORLD, "\n\n*** End of run ***\n\n"); CHKERRQ(ierr);

	ierr = PetscFinalize(); CHKERRQ(ierr);

	return 0;
} // main
