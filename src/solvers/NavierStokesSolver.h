#if !defined(NAVIER_STOKES_SOLVER_H)
#define NAVIER_STOKES_SOLVER_H

#include "FlowDescription.h"
#include "CartesianMesh.h"
#include "SimulationParameters.h"
//#include "FractionalStepMethod.h"
#include <petscdmda.h>
#include <petscksp.h>

template <PetscInt dim>
class NavierStokesSolver
{
protected:
	// classes
	FlowDescription      *flowDesc;
	SimulationParameters *simParams;
	CartesianMesh        *mesh;
	//FractionalStepMethod FSM;
	
	PetscInt timeStep,
	         iteratonCount1,
	         iterationCount2;

	std::vector<PetscReal> dxU, dyU, dzU,
	                       dxV, dyV, dzV,
	                       dxW, dyW, dzW;
	
	DM  pda,
	    uda,
	    vda,
	    wda,
	    qPack,
	    phiPack;
	
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
	Vec bc1, rhs1, rhs2;
	Vec q, qStar, phi;
	KSP ksp1, ksp2;

	void createDMs();
	void createVecs();
	void createKSPs();
	void initialiseMeshSpacings();
	void initialiseFluxes();
	void createLocalToGlobalMappingsFluxes();
	void createLocalToGlobalMappingsPhi();
	void updateBoundaryGhosts();
	void generateDiagonalMatrices();
	void generateA();
	void calculateExplicitTerms();
	void generateBC1();
	void generateRHS1();
	void generateBNQ();
	void solveIntermediateVelocity();
	void projectionStep();
	
public:
	void initialise();
	void finalise();
	void stepTime();
	void writeData(std::string caseFolder);
	bool finished();
	
	/**
	* @brief Give the name of the current solver 
	* @return String that describes the type of solver
	*/
	virtual std::string name()
	{
		return "Navier-Stokes";
	}
	
	NavierStokesSolver(FlowDescription *FD, SimulationParameters *SP, CartesianMesh *CM)
	{
		timeStep = 0;
		// classes
		flowDesc  = FD;
		simParams = SP;
		mesh      = CM;
		// DMs
		pda = PETSC_NULL;
		uda = PETSC_NULL;
		vda = PETSC_NULL;
		wda = PETSC_NULL;
		qPack   = PETSC_NULL;
		phiPack = PETSC_NULL;
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
	}
};

#endif