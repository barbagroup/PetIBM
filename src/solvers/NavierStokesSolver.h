/***************************************************************************//**
* \file
* \brief Header to define class NavierStokesSolver.
*/

#if !defined(NAVIER_STOKES_SOLVER_H)
#define NAVIER_STOKES_SOLVER_H

#include "FlowDescription.h"
#include "CartesianMesh.h"
#include "SimulationParameters.h"
#include <petscdmda.h>
#include <petscksp.h>
#include <fstream>

/***************************************************************************//**
* \brief Solve the incompressible Navier-Stokes equations in a rectangular or
*        cuboidal domain.
*/
template <PetscInt dim>
class NavierStokesSolver
{
protected:
	std::string caseFolder;

	FlowDescription      *flowDesc;
	SimulationParameters *simParams;
	CartesianMesh        *mesh;
	
	PetscInt timeStep,
	         iteratonCount1,
	         iterationCount2;

	std::vector<PetscReal> dxU, dyU, dzU,
	                       dxV, dyV, dzV,
	                       dxW, dyW, dzW;

	std::ofstream iterationsFile;
	
	DM  pda,
	    uda,
	    vda,
	    wda,
	    qPack,
	    lambdaPack;
	
	Vec qxLocal,
	    qyLocal,
	    qzLocal;

	Vec uMapping,
	    vMapping,
	    wMapping,
	    pMapping;

	Vec H, rn;
	Vec RInv, MHat;

	Mat A;
	Mat QT, BNQ;
	Mat QTBNQ;
	Vec BN;
	Vec bc1, rhs1, r2, rhs2, temp;
	Vec q, qStar, lambda;
	KSP ksp1, ksp2;
	PC  pc2;

	PetscLogStage stageInitialize,
	              stageSolveIntermediateVelocity,
	              stageSolvePoissonSystem,
	              stageProjectionStep;

	/**
	* \brief Initialize data common to NavierStokesSolver and derived classes
	*/
	PetscErrorCode initializeCommon();

	/**
	* \brief Create the DMDA structures for the flow variables
	*/
	virtual PetscErrorCode createDMs();

	/**
	* \brief Create the vectors used to store the flow variables
	*/
	virtual PetscErrorCode createVecs();

	/**
	* \brief Set up the Kyrlov solvers used to solve the linear systems
	*/
	PetscErrorCode createKSPs();

	/**
	* \brief Initialize the spaces between adjacent velocity nodes
	*/
	void initializeMeshSpacings();

	/**
	* \brief Populate flux vectors with the initial conditions
	*/
	PetscErrorCode initializeFluxes();

	/**
	* \brief Read the fluxes from previously saved data
	*/
	PetscErrorCode readFluxes();

	/**
	* \brief Initialize lamba with previously saved data
	*/
	virtual PetscErrorCode initializeLambda();

	/**
	* \brief Create the mappings from the local flux variables to the global 
	*        flux vector
	*/
	PetscErrorCode createLocalToGlobalMappingsFluxes();

	/**
	* \brief Create the mapping from the local pressure variables to the global 
	*        lambda vector
	*/
	PetscErrorCode createLocalToGlobalMappingsLambda();

	/**
	* \brief Update the values in the ghost nodes on the domain boundary
	*/
	PetscErrorCode updateBoundaryGhosts();

	/**
	* \brief Generate the diagonal matrices M and Rinv
	*/
	PetscErrorCode generateDiagonalMatrices();

	/**
	* \brief Count the number of non-zeros in the diagonal and off-diagonal 
	*        portions of the parallel matrices
	*/
	void countNumNonZeros(PetscInt *cols, size_t numCols, PetscInt rowStart, PetscInt rowEnd, PetscInt &d_nnz, PetscInt &o_nnz);

	/**
	* \brief Generate the matrix \f$ A \f$
	*/
	PetscErrorCode generateA();

	/**
	* \brief Calculate the explicit convection and diffusion terms
	*/
	PetscErrorCode calculateExplicitTerms();

	/**
	* \brief Assemble the vector arising from the boundary conditons in the RHS 
	*        of the intermediate velocity solve
	*/
	PetscErrorCode generateBC1();

	/**
	* \brief Calculate the RHS of the intermediate velocity solve
	*/
	PetscErrorCode generateRHS1();

	/**
	* \brief Assemble the vector arising from the boundary conditions in the 
	*        RHS of the pressure-force solve
	*/
	virtual PetscErrorCode generateR2();

	/**
	* \brief Generate the RHS of the pressure-force solve
	*/
	PetscErrorCode generateRHS2();

	/**
	* \brief Assemble the matrix \f$ B^N Q \f$
	*/
	virtual PetscErrorCode generateBNQ();

	/**
	* \brief Calculate the matrix \f$ Q^T B^N Q \f$
	*/
	PetscErrorCode generateQTBNQ();

	/**
	* \brief Calculate and specify to the Krylov solver the null space of the 
	*        LHS matrix in the pressure-force solve
	*/
	virtual PetscErrorCode setNullSpace();

	/**
	* \brief Solve for the intermediate velocity fluxes \f$ q^* \f$
	*/
	PetscErrorCode solveIntermediateVelocity();

	/**
	* \brief Solve for the pressure and body forces
	*/
	PetscErrorCode solvePoissonSystem();

	/**
	* \brief Project the pressure and forces on to the velocity field to obtain 
	*        the velocity at the next time step
	*/
	PetscErrorCode projectionStep();

	/**
	* \brief Write the velocity fluxes to files
	*/
	PetscErrorCode writeFluxes();

	/**
	* \brief Write the pressure field to file
	*/
	virtual PetscErrorCode writeLambda();
	
public:
	/**
	* \brief Initial set-up of the system
	*/
	virtual PetscErrorCode initialize();

	/**
	* \brief Clean up at the end of the simulation
	*/
	virtual PetscErrorCode finalize();

	/**
	* \brief Move the simulation forward one time step
	*/
	PetscErrorCode stepTime();

	/**
	* Write the flow variables to files
	*/
	virtual PetscErrorCode writeData();

	/**
	* \brief Write the simulation parameters to file
	*/
	PetscErrorCode writeSimulationInfo();

	/**
	* \brief Write the grid coordinates to file
	*/
	PetscErrorCode writeGrid();

	/**
	* \brief Specify if the data needs to be saved in the current time step
	*/
	PetscBool savePoint();

	/**
	* \brief Specify if the simulation is completed.
	*/
	PetscBool finished();
	
	/**
	* \brief Give the name of the current solver
	* \return String that describes the type of solver
	*/
	virtual std::string name()
	{
		return "Navier-Stokes";
	}
	
	NavierStokesSolver(std::string folder, FlowDescription *FD, SimulationParameters *SP, CartesianMesh *CM)
	{
		// classes
		caseFolder= folder;
		flowDesc  = FD;
		simParams = SP;
		mesh      = CM;
		timeStep  = simParams->startStep;
		// DMs
		pda = PETSC_NULL;
		uda = PETSC_NULL;
		vda = PETSC_NULL;
		wda = PETSC_NULL;
		qPack   = PETSC_NULL;
		lambdaPack = PETSC_NULL;
		// Vecs
		qxLocal  = PETSC_NULL;
		qyLocal  = PETSC_NULL;
		qzLocal  = PETSC_NULL;
		q        = PETSC_NULL;
		qStar    = PETSC_NULL;
		H        = PETSC_NULL;
		rn       = PETSC_NULL;
		bc1      = PETSC_NULL;
		rhs1     = PETSC_NULL;
		r2       = PETSC_NULL;
		rhs2     = PETSC_NULL;
		temp     = PETSC_NULL;
		RInv     = PETSC_NULL;
		MHat	 = PETSC_NULL;
		BN       = PETSC_NULL;
		pMapping = PETSC_NULL;
		uMapping = PETSC_NULL;
		vMapping = PETSC_NULL;
		wMapping = PETSC_NULL;
		// Mats
		A       = PETSC_NULL;
		QT      = PETSC_NULL;
		BNQ     = PETSC_NULL;
		QTBNQ   = PETSC_NULL;
		//KSPs
		ksp1 = PETSC_NULL;
		ksp2 = PETSC_NULL;
		// PCs
		pc2 = PETSC_NULL;
		// PetscLogStages
		PetscLogStageRegister("initialize", &stageInitialize);
		PetscLogStageRegister("solveIntVel", &stageSolveIntermediateVelocity);
		PetscLogStageRegister("solvePoissSys", &stageSolvePoissonSystem);
		PetscLogStageRegister("projectionStep", &stageProjectionStep);
	}
};

#endif
