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
  Mat E,          ///< interpolation operator
      ET,         ///< spreading operator
      EBNET,      ///< matrix of the system for the Lagrangian forces
      G;          ///< gradient operator
  Vec fTilde,     ///< vector for the Lagrangian forces
      rhsf,       ///< right-hand side of the system for the Lagrangian forces
      dfTilde;    ///< delta (variation) of the Lagrangian forces
  Vec dlambda;    ///< delta (variation) of the pressure field
  Vec tmp;        ///< a temporary vector
  Vec rhs1_n;     ///< RHS of velocity system without pressure gradient and momentum forcing
  Solver *forces; ///< solver for the Lagrangian forces

  PetscReal bodyForces[dim]; ///< array with the force in each direction acting on the body
  std::ofstream forcesFile;  ///< file in which to write the forces

  PetscInt algorithm;      ///< algo index for order of decoupling
  PetscInt forceEstimator; ///< scheme index to estimate momentum forcing
  PetscReal atol,          ///< absolute tolerance for sub-iterative process
            rtol;          ///< relative tolerance for sub-iterative process
  PetscInt maxIters;       ///< maximum number of iterations for sub-iterative process
  PetscBool printStats;    ///< prints L2 norm and relative L2 norm when using sub-iterative process

  PetscLogStage stageRHSForceSystem,
                stageSolveForceSystem,
                stageIntegrateForces;

  PetscErrorCode initializeBodies();
  PetscErrorCode getNumLagPoints(PetscInt &n);
  PetscErrorCode getNumLagPointsOnProcess(std::vector<PetscInt> &numOnProcess);
  PetscErrorCode registerLagPointsOnProcess();
  PetscErrorCode createDMs();
  PetscErrorCode createVecs();
  PetscErrorCode createGlobalMappingBodies();
  PetscErrorCode generateET();
  PetscErrorCode updateRHSVelocity();
  PetscErrorCode solvePoissonSystem(Vec &p);
  PetscErrorCode createForceSolver();
  PetscErrorCode updateFlux(Vec f);
  PetscErrorCode assembleRHSForce(Vec q);
  PetscErrorCode solveForceSystem(Vec &f);
  PetscErrorCode projectionStep(Vec p);

  PetscErrorCode calculateForces();
  PetscErrorCode calculateForces2();
  PetscErrorCode writeForces();
  PetscErrorCode writeIterationCounts();
  PetscErrorCode writeLagrangianForces(std::string directory);

  PetscErrorCode scatterGlobalToLocal();

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
