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

	std::vector<PetscReal> x, y, z;
	std::vector<PetscInt>  globalIndexMapping;
	std::vector<PetscInt>  startGlobalIndices;
	std::vector<PetscInt>  numBoundaryPointsOnProcess;
	std::vector<PetscInt>  numPhiOnProcess;
	std::vector< std::vector<PetscInt> > boundaryPointIndices;
	
	PetscErrorCode createDMs();
	PetscErrorCode generateBNQ();
	PetscErrorCode generateET();
	PetscErrorCode generateR2();
	PetscErrorCode createGlobalMappingBodies();

	PetscReal dhRoma(PetscReal x, PetscReal h);
	PetscReal delta(PetscReal x, PetscReal y, PetscReal h);
	PetscReal delta(PetscReal x, PetscReal y, PetscReal z, PetscReal h);

public:
	PetscErrorCode initialise();
	void initialiseBodies();

	TairaColoniusSolver(std::string folder, FlowDescription *FD, SimulationParameters *SP, CartesianMesh *CM) : NavierStokesSolver<dim>::NavierStokesSolver(folder, FD, SP, CM)
	{
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