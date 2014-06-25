#if !defined(TAIRA_COLONIUS_SOLVER_H)
#define TAIRA_COLONIUS_SOLVER_H

#include "NavierStokesSolver.h"

template <PetscInt dim>
class TairaColoniusSolver : public NavierStokesSolver<dim>
{
protected:
	PetscInt  startGlobalIndex;
	DM        bda;
	Mat       ET;
	PetscReal force[3];
	Vec       temp;

	std::ofstream forcesFile;

	std::vector<PetscReal> x, y, z;
	std::vector<PetscInt>  I, J, K;
	std::vector<PetscInt>  globalIndexMapping;
	std::vector<PetscInt>  numBoundaryPointsOnProcess;
	std::vector<PetscInt>  numPhiOnProcess;
	std::vector< std::vector<PetscInt> > boundaryPointIndices;
	
	PetscErrorCode initializeLambda();
	PetscErrorCode initializeBodies();
	PetscErrorCode generateBodyInfo();
	PetscErrorCode calculateCellIndices();
	PetscErrorCode createDMs();
	PetscErrorCode generateBNQ();
	PetscErrorCode generateR2();
	PetscErrorCode createGlobalMappingBodies();
	PetscErrorCode calculateForce();
	PetscErrorCode writeForces();
	PetscErrorCode writeLambda();

	PetscReal dhRoma(PetscReal x, PetscReal h);
	PetscReal delta(PetscReal x, PetscReal y, PetscReal h);
	PetscReal delta(PetscReal x, PetscReal y, PetscReal z, PetscReal h);
	PetscBool isInfluenced(PetscReal xGrid, PetscReal yGrid, PetscReal xBody, PetscReal yBody, PetscReal radius, PetscReal *delta);
	PetscBool isInfluenced(PetscReal xGrid, PetscReal yGrid, PetscReal zGrid, PetscReal xBody, PetscReal yBody, PetscReal zBody, PetscReal radius, PetscReal *delta);

public:
	PetscErrorCode initialize();
	PetscErrorCode finalize();
	PetscErrorCode writeData();

	TairaColoniusSolver(std::string folder, FlowDescription *FD, SimulationParameters *SP, CartesianMesh *CM) : NavierStokesSolver<dim>::NavierStokesSolver(folder, FD, SP, CM)
	{
		bda = PETSC_NULL;
		ET  = PETSC_NULL;
		temp= PETSC_NULL;
	}
	
	/**
	* @brief Give the name of the current solver 
	* @return String that describes the type of solver
	*/
	virtual std::string name()
	{
		return "Taira and Colonius";
	}
};

#endif