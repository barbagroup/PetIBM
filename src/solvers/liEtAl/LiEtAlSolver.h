/*! Definition of the class `LiEtAlSolver`.
 * \file LiEtAlSolver.h
 */


#if !defined(LI_ET_AL_SOLVER_H)
#define LI_ET_AL_SOLVER_H

#include "navierStokes/NavierStokesSolver.h"

#include "Body.h"


/*!
 * \class LiEtAlSolver
 * \brief Solves the Navier-Stokes equations 
 *        with immersed boundary projection method (Li et al. , 2016).
 */
template <PetscInt dim>
class LiEtAlSolver : public NavierStokesSolver<dim>
{
public:
  PetscInt numBodies;
  std::vector<Body<dim> > bodies;

  DM bda;
  Mat E, ET, EBNET;
  Vec fTilde, rhsf, temp, temp2, qStar2;
  Solver *forces;

  std::ofstream forcesFile;

  PetscLogStage stageRHSForceSystem,
                stageSolveForceSystem;

  PetscErrorCode initializeBodies();
  PetscErrorCode getNumLagPoints(PetscInt &n);
  PetscErrorCode getNumLagPointsOnProcess(std::vector<PetscInt> &numOnProcess);
  PetscErrorCode registerLagPointsOnProcess();
  PetscErrorCode createDMs();
  PetscErrorCode createVecs();
  PetscErrorCode createGlobalMappingBodies();
  PetscErrorCode generateET();
  PetscErrorCode createForceSolver();
  PetscErrorCode updateFlux();
  PetscErrorCode assembleRHSForce();
  PetscErrorCode solveForceSystem();

  PetscErrorCode calculateForces();
  PetscErrorCode writeForces();
  PetscErrorCode writeIterationCounts();
  PetscErrorCode writeLagrangianForces(std::string directory);

public:
  // constructors
  LiEtAlSolver(){ };
  LiEtAlSolver(CartesianMesh *cartesianMesh, 
               FlowDescription<dim> *flowDescription, 
               SimulationParameters *simulationParameters);
  // destructor
  ~LiEtAlSolver(){ };
  PetscErrorCode initialize();
  PetscErrorCode stepTime();
  PetscErrorCode finalize();
  PetscErrorCode writeData();

}; // LiEtAlSolver

#endif
