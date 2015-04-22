/***************************************************************************//**
 * \file NavierStokesSolver.h
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Definition of the class \c NavierStokesSolver.
 */


#if !defined(NAVIER_STOKES_SOLVER_H)
#define NAVIER_STOKES_SOLVER_H

#include "FlowDescription.h"
#include "CartesianMesh.h"
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
  std::string caseFolder;

  FlowDescription      *flowDesc;
  SimulationParameters *simParams;
  CartesianMesh        *mesh;
  
  PetscInt timeStep,
           iteratonCount1,
           iterationCount2;

  std::vector<PetscReal> dxU, dyU, dzU,
                         dxV, dyV, dzV,
                         dxW, dyW, dzW;

  std::ofstream iterationsFile;
  
  DM  pda,
      uda,
      vda,
      wda,
      qPack,
      lambdaPack;
  
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
  Vec bc1, rhs1, r2, rhs2, temp;
  Vec q, qStar, lambda;
  KSP ksp1, ksp2;
  PC  pc2;

  PetscLogStage stageInitialize,
                stageSolveIntermediateVelocity,
                stageSolvePoissonSystem,
                stageProjectionStep;

  // initialize data common to NavierStokesSolver and derived classes
  PetscErrorCode initializeCommon();

  // create DMDA structures for flow variables
  virtual PetscErrorCode createDMs();

  // create vectors used to store flow variables
  virtual PetscErrorCode createVecs();

  // set up Krylov solvers used to solve linear systems
  PetscErrorCode createKSPs();

  // initialize spaces between adjacent velocity nodes
  void initializeMeshSpacings();

  // populate flux vectors with initial conditions
  virtual PetscErrorCode initializeFluxes();

  // read fluxes from previously saved data
  PetscErrorCode readFluxes();

  // initialize lambda vector with previously saved data
  virtual PetscErrorCode initializeLambda();

  // create mapping from local flux vectors to global flux vectors
  PetscErrorCode createLocalToGlobalMappingsFluxes();

  // create mapping from local pressure variable to global lambda vector
  PetscErrorCode createLocalToGlobalMappingsLambda();

  // update values in ghost nodes at the domain boundaries
  PetscErrorCode updateBoundaryGhosts();

  // generate diagonal matrices M and Rinv
  PetscErrorCode generateDiagonalMatrices();

  // count number of non-zeros in the diagonal and off-diagonal portions of the parallel matrices
  void countNumNonZeros(PetscInt *cols, size_t numCols, PetscInt rowStart, PetscInt rowEnd, PetscInt &d_nnz, PetscInt &o_nnz);

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

  // solve system for intermediate velocity fluxes \f$ q^* \f$
  PetscErrorCode solveIntermediateVelocity();

  // solver Poisson system for pressure and body forces
  PetscErrorCode solvePoissonSystem();

  // project velocity onto divergence-free field with satisfaction of the no-splip condition
  PetscErrorCode projectionStep();

  // write fluxes into files
  PetscErrorCode writeFluxes();

  // write pressure filed into file
  virtual PetscErrorCode writeLambda();
  
public:
  // initial set-up of the system
  virtual PetscErrorCode initialize();

  // clean up at the end of the simulation
  virtual PetscErrorCode finalize();

  // advance in time
  PetscErrorCode stepTime();

  // write flow variables into files
  virtual PetscErrorCode writeData();

  // write simulation parameters into file
  PetscErrorCode printSimulationInfo();

  // write grid coordinates into file
  PetscErrorCode writeGrid();

  // specify if data needs to be saved at current time-step
  PetscBool savePoint();

  // evaluate if the simulation is completed
  PetscBool finished();
  
  // name of the solver
  virtual std::string name()
  {
    return "Navier-Stokes";
  }
  
  NavierStokesSolver(std::string folder, FlowDescription *FD, SimulationParameters *SP, CartesianMesh *CM)
  {
    // classes
    caseFolder= folder;
    flowDesc  = FD;
    simParams = SP;
    mesh      = CM;
    timeStep  = simParams->startStep;
    // DMs
    pda = PETSC_NULL;
    uda = PETSC_NULL;
    vda = PETSC_NULL;
    wda = PETSC_NULL;
    qPack   = PETSC_NULL;
    lambdaPack = PETSC_NULL;
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
    r2       = PETSC_NULL;
    rhs2     = PETSC_NULL;
    temp     = PETSC_NULL;
    RInv     = PETSC_NULL;
    MHat   = PETSC_NULL;
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
    // PCs
    pc2 = PETSC_NULL;
    // PetscLogStages
    PetscLogStageRegister("initialize", &stageInitialize);
    PetscLogStageRegister("solveIntVel", &stageSolveIntermediateVelocity);
    PetscLogStageRegister("solvePoissSys", &stageSolvePoissonSystem);
    PetscLogStageRegister("projectionStep", &stageProjectionStep);
  }
};

#endif