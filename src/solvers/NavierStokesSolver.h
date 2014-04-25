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
	Vec H, rn;
	Mat Rinv, M;

	void fluxVecsCreate();
	void fluxVecsInitialise();
	void updateBoundaryGhosts();
	void calculateExplicitTerms();
	
public:
	void initialise();
	void stepTime();
	void writeData();
	bool finished();
	
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
	
	NavierStokesSolver()
	{
		timeStep = 0;
		flowDesc  = NULL;
		simParams = NULL;
		mesh      = NULL;
		uda  = PETSC_NULL;
		vda  = PETSC_NULL;
		wda  = PETSC_NULL;
		pack = PETSC_NULL;
		qxLocal = PETSC_NULL;
		qyLocal = PETSC_NULL;
		qzLocal = PETSC_NULL;
		H  = PETSC_NULL;
		rn = PETSC_NULL;
		Rinv = PETSC_NULL;
		M    = PETSC_NULL;
	}
	NavierStokesSolver(FlowDescription &FD, SimulationParameters &SP, CartesianMesh &CM)
	{
		timeStep = 0;
		flowDesc  = &FD;
		simParams = &SP;
		mesh      = &CM;
		uda  = PETSC_NULL;
		vda  = PETSC_NULL;
		wda  = PETSC_NULL;
		pack = PETSC_NULL;
		qxLocal = PETSC_NULL;
		qyLocal = PETSC_NULL;
		qzLocal = PETSC_NULL;
		H  = PETSC_NULL;
		rn = PETSC_NULL;
		Rinv = PETSC_NULL;
		M    = PETSC_NULL;
	}
	virtual ~NavierStokesSolver()
	{
		PetscErrorCode ierr;
		if(!H)
		{
			ierr = VecDestroy(&H); CHKERRV(ierr);
		}
		if(!rn)
		{
			ierr = VecDestroy(&rn); CHKERRV(ierr);
		}
		if(!qxLocal)
		{
			ierr = VecDestroy(&qxLocal); CHKERRV(ierr);
		}
		if(!qyLocal)
		{
			ierr = VecDestroy(&qyLocal); CHKERRV(ierr);
		}
		if(!qzLocal)
		{
			ierr = VecDestroy(&qzLocal); CHKERRV(ierr);
		}
		if(!uda)
		{
			ierr = DMDestroy(&uda); CHKERRV(ierr);
		}
		if(!vda)
		{
			ierr = DMDestroy(&vda); CHKERRV(ierr);
		}
		if(!wda)
		{
			ierr = DMDestroy(&wda); CHKERRV(ierr);
		}
		if(!pack)
		{
			ierr = DMDestroy(&pack); CHKERRV(ierr);
		}
		if(flowDesc!=NULL)
			delete flowDesc;
		if(simParams!=NULL)
			delete simParams;
		if(mesh!=NULL)
			delete mesh;
	}
};

#endif
