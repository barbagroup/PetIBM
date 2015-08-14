/***************************************************************************//**
 * \file NavierStokesSolver.h
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Definition of the class \c NavierStokesSolver.
 */


#if !defined(NAVIER_STOKES_SOLVER_H)
#define NAVIER_STOKES_SOLVER_H

#include "CartesianMesh.h"
#include "FlowDescription.h"
#include "SimulationParameters.h"

#include <fstream>

#include <petscdmda.h>
#include <petscksp.h>


/**
 * \brief Solve the incompressible Navier-Stokes equations in a rectangular or
 *        cuboidal domain.
 */
template <PetscInt dim>
class NavierStokesSolver
{
public:
  DM qPack,
     lambdaPack;
  DM pda,
     uda,
     vda,
     wda;

  Vec q, qStar, lambda;

  Vec qxLocal, qyLocal, qzLocal;
  Vec pMapping, uMapping, vMapping, wMapping;

  KSP ksp1, ksp2;
  PC pc1, pc2;

  Mat A,
      QT,
      BNQ,
      QTBNQ;

  Vec bc1,
      rhs1,
      H,
      rn,
      r2,
      rhs2,
      temp;

  Vec BN,
      RInv,
      MHat;

  CartesianMesh *mesh;
  FlowDescription<dim> *flow;
  SimulationParameters *parameters;
  
  PetscInt timeStep;

  std::vector<PetscReal> dxU, dyU, dzU,
                         dxV, dyV, dzV,
                         dxW, dyW, dzW;

  std::ofstream iterationCountsFile;

  PetscLogStage stageInitialize,
                stageRHSVelocitySystem,
                stageSolveVelocitySystem,
                stageRHSPoissonSystem,
                stageSolvePoissonSystem,
                stageProjectionStep;

  // initialize data common to NavierStokesSolver and derived classes
  PetscErrorCode initializeCommon();
  // create DMDA structures for flow variables
  virtual PetscErrorCode createDMs();
  // create vectors used to store flow variables
  virtual PetscErrorCode createVecs();
  // create mapping from local flux vectors to global flux vectors
  PetscErrorCode createLocalToGlobalMappingsFluxes();
  // create mapping from local pressure variable to global lambda vector
  PetscErrorCode createLocalToGlobalMappingsLambda();
  // set up Krylov solvers used to solve linear systems
  PetscErrorCode createSolvers();
  // set up Krylov solver for velocity system
  PetscErrorCode createVelocitySolver();
  // set up Krylov solver for Poisson system
  PetscErrorCode createPoissonSolver();
  // initialize spaces between adjacent velocity nodes
  void initializeMeshSpacings();
  // populate flux vectors with initial conditions
  PetscErrorCode initializeFluxes();
  PetscErrorCode addInitialPerturbation();
  // initialize lambda vector with previously saved data
  PetscErrorCode initializeLambda();
  // generate diagonal matrices M and Rinv
  PetscErrorCode generateDiagonalMatrices();
  // count number of non-zeros in the diagonal and off-diagonal portions of the parallel matrices
  void countNumNonZeros(PetscInt *cols, size_t numCols, PetscInt rowStart, PetscInt rowEnd, 
                        PetscInt &d_nnz, PetscInt &o_nnz);
  // generate the matrix A
  PetscErrorCode generateA();
  // calculate explicit convective and diffusive terms
  PetscErrorCode calculateExplicitTerms();
  // assemble velocity boundary conditions vector
  PetscErrorCode generateBC1();
  // assemble RHS of intermediate velocity system
  PetscErrorCode generateRHS1();
  // assemble boundary conditions vector for the pressure-force system
  virtual PetscErrorCode generateR2();
  // assemble RHS of pressure-force system
  PetscErrorCode generateRHS2();
  // compute matrix \f$ B^N Q \f$
  virtual PetscErrorCode generateBNQ();
  // compute matrix \f$ Q^T B^N Q \f$
  PetscErrorCode generateQTBNQ();
  // calculate and specify to the Krylov solver the null-space of the LHS matrix
  // in the pressure-force system
  virtual PetscErrorCode setNullSpace();
  // update values in ghost nodes at the domain boundaries
  PetscErrorCode updateBoundaryGhosts();

  // advance in time
  PetscErrorCode stepTime();
  // solve system for intermediate velocity fluxes \f$ q^* \f$
  PetscErrorCode solveIntermediateVelocity();
  // solver Poisson system for pressure and body forces
  PetscErrorCode solvePoissonSystem();
  // project velocity onto divergence-free field with satisfaction of the no-splip condition
  PetscErrorCode projectionStep();

  // print info about simulation
  PetscErrorCode printInfo();
  // read fluxes from files
  PetscErrorCode readFluxes();
  // read pressure field from file
  virtual PetscErrorCode readLambda();
  // write grid coordinates into file
  PetscErrorCode writeGrid();
  // write fluxes into files
  PetscErrorCode writeFluxes();
  // write pressure field into file
  virtual PetscErrorCode writeLambda();
  // write KSP iteration counts into file
  PetscErrorCode writeIterationCounts();
  
public:
  // constructor
  NavierStokesSolver(CartesianMesh *cartesianMesh, 
                     FlowDescription<dim> *flowDescription, 
                     SimulationParameters *simulationParameters);
  // destructor
  ~NavierStokesSolver();

  // initialize systems
  PetscErrorCode initialize();
  // clean up at end of simulation
  PetscErrorCode finalize();

  // specify if data needs to be saved at current time-step
  PetscBool savePoint();
  // evaluate if the simulation is completed
  PetscBool finished();

  // write numerical solution into respective files
  virtual PetscErrorCode writeData();

}; // NavierStokesSolver

#endif