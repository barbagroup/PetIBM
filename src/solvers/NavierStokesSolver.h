#if !defined(NAVIER_STOKES_SOLVER_H)
#define NAVIER_STOKES_SOLVER_H

#include "FlowDescription.h"
#include "CartesianMesh.h"
#include "SimulationParameters.h"
#include "FractionalStepMethod.h"
#include <petscdmda.h>

template <PetscInt dim>
class NavierStokesSolver
{
protected:
	// classes
	FlowDescription      *flowDesc;
	SimulationParameters *simParams;
	CartesianMesh        *mesh;
	FractionalStepMethod FSM;
	
	size_t timeStep, iteratonCount1, iterationCount2;
	
	DM  uda, vda, wda, pack;
	Vec qxLocal, qyLocal, qzLocal;
	Mat Rinv, M;

	void fluxVecsCreate();
	void fluxVecsInitialise();
	void updateBoundaryGhosts();
	
public:
	void initialise();
	
	// Factory methods are static (not entirely sure why)
	static NavierStokesSolver<dim>* createSolver(FlowDescription &FD, SimulationParameters &SP, CartesianMesh &CM);
	
	/**
	* @brief Give the name of the current solver 
	* @return String that describes the type of solver
	*/
	virtual std::string name()
	{
		return "Navier-Stokes";
	}
};

#endif
