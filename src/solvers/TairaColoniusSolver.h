#if !defined(TAIRA_COLONIUS_SOLVER_H)
#define TAIRA_COLONIUS_SOLVER_H

#include "NavierStokesSolver.h"

template <PetscInt dim>
class TairaColoniusSolver : public NavierStokesSolver<dim>
{
protected:
	DM        bda;
	PetscInt  startGlobalIndex;
	PetscReal h;

	std::vector<PetscReal> x, y, z;
	std::vector<PetscReal> ds;
	std::vector<PetscInt>  bodyGlobalIndices;
	std::vector<PetscInt>  startGlobalIndices;
	std::vector<PetscInt>  numBoundaryPointsOnProcess;
	std::vector< std::vector<PetscInt> > boundaryPointIndices;
	
	void generateBNQ();
	void generateR2();

public:
	void initialise();
	void initialiseBodies();
	void finalise();

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