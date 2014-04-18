#if !defined(NAVIER_STOKES_SOLVER_H)
#define NAVIER_STOKES_SOLVER_H

#include "CartesianMesh.h"
#include "SimulationParameters.h"

template <int dim>
class NavierStokesSolver
{
protected:
	// classes
	SimulationParameters *simParams;
	CartesianMesh<dim>   *mesh;
	
public:	
	// Factory methods are static (not entirely sure why)
	static NavierStokesSolver<dim>* createSolver(SimulationParameters &SP, CartesianMesh<dim> &CM);
	
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
