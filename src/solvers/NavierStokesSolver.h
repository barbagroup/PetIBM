#if !defined(NAVIER_STOKES_SOLVER_H)
#define NAVIER_STOKES_SOLVER_H

#include "FlowDescription.h"
#include "CartesianMesh.h"
#include "SimulationParameters.h"
#include <petscdmda.h>
#include <petscksp.h>
#include <fstream>

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

	PetscErrorCode initializeCommon();
	virtual PetscErrorCode createDMs();
	virtual PetscErrorCode createVecs();
	PetscErrorCode createKSPs();
	void initializeMeshSpacings();
	PetscErrorCode initializeFluxes();
	PetscErrorCode readFluxes(Vec qxGlobal, Vec qyGlobal, Vec qzGlobal=PETSC_NULL);
	virtual PetscErrorCode initializeLambda();
	PetscErrorCode createLocalToGlobalMappingsFluxes();
	PetscErrorCode createLocalToGlobalMappingsLambda();
	PetscErrorCode updateBoundaryGhosts();
	PetscErrorCode generateDiagonalMatrices();
	void countNumNonZeros(PetscInt *cols, size_t numCols, PetscInt rowStart, PetscInt rowEnd, PetscInt &d_nnz, PetscInt &o_nnz);
	PetscErrorCode generateA();
	PetscErrorCode calculateExplicitTerms();
	PetscErrorCode generateBC1();
	PetscErrorCode generateRHS1();
	virtual PetscErrorCode generateR2();
	PetscErrorCode generateRHS2();
	virtual PetscErrorCode generateBNQ();
	PetscErrorCode generateQTBNQ();
	virtual PetscErrorCode setNullSpace();
	PetscErrorCode solveIntermediateVelocity();
	PetscErrorCode solvePoissonSystem();
	PetscErrorCode projectionStep();
	PetscErrorCode writeFluxes();
	virtual PetscErrorCode writeLambda();
	
public:
	virtual PetscErrorCode initialize();
	virtual PetscErrorCode finalize();
	PetscErrorCode stepTime();
	virtual PetscErrorCode writeData();
	PetscErrorCode writeSimulationInfo();
	PetscErrorCode writeGrid();
	PetscBool savePoint();
	PetscBool finished();
	
	/**
	* @brief Give the name of the current solver 
	* @return String that describes the type of solver
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
