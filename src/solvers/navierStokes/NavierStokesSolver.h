/***************************************************************************//**
 * \file NavierStokesSolver.h
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Definition of the class `NavierStokesSolver`.
 */


#if !defined(NAVIER_STOKES_SOLVER_H)
#define NAVIER_STOKES_SOLVER_H

#include "CartesianMesh.h"
#include "FlowDescription.h"
#include "SimulationParameters.h"
#include "solvers/solver.h"

#include <fstream>
#include <memory>

#include <petscdmda.h>
#include <petscksp.h>


/**
 * \class NavierStokesSolver
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

  Solver *velocity, *poisson;

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
  // create the solvers for the linear systems
  PetscErrorCode createSolvers();
  // populate flux vectors with initial conditions
  virtual PetscErrorCode initializeFluxes();
  // add initial perturbation to velocity field
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
  
  // compute matrix \f$ B^N Q \f$
  virtual PetscErrorCode generateBNQ();
  // generate gradient operator
  PetscErrorCode generateGradient(Mat *G);
  // compute matrix \f$ Q^T B^N Q \f$
  PetscErrorCode generateQTBNQ();
  // calculate and specify to the Krylov solver the null-space of the LHS matrix
  // in the pressure-force system
  virtual PetscErrorCode setNullSpace();

  // assemble RHS of velocity system
  PetscErrorCode assembleRHSVelocity();
  // calculate explicit convective and diffusive terms
  PetscErrorCode calculateExplicitTerms();
  // update values in ghost nodes at the domain boundaries
  PetscErrorCode updateBoundaryGhosts();
  // assemble velocity boundary conditions vector
  PetscErrorCode generateBC1();

  // aasemble RHS of Poisson system
  PetscErrorCode assembleRHSPoisson();
  // assemble boundary conditions vector for the pressure-force system
  virtual PetscErrorCode generateR2();

  // advance in time
  virtual PetscErrorCode stepTime();
  // solve system for intermediate velocity fluxes \f$ q^* \f$
  PetscErrorCode solveIntermediateVelocity();
  // solver Poisson system for pressure and body forces
  virtual PetscErrorCode solvePoissonSystem();
  // project velocity onto divergence-free field with satisfaction of the no-splip condition
  virtual PetscErrorCode projectionStep();

  // read fluxes from files
  PetscErrorCode readFluxes(std::string directory);
  // read velocity components from files
  PetscErrorCode readVelocities(std::string directory);
  // read convective terms from files
  PetscErrorCode readConvectiveTerms(std::string directory);
  // read pressure field from file
  virtual PetscErrorCode readLambda(std::string directory);
#ifdef PETSC_HAVE_HDF5
  // write grid stations of the different field variables into HDF5 files
  PetscErrorCode writeGrids();
#endif
  // write fluxes into files
  PetscErrorCode writeFluxes(std::string directory);
  // write velocity components into files
  PetscErrorCode writeVelocities(std::string directory);
  // write convective terms into files
  PetscErrorCode writeConvectiveTerms(std::string directory);
  // write pressure field into file
  virtual PetscErrorCode writeLambda(std::string directory);
  // write KSP iteration counts into file
  virtual PetscErrorCode writeIterationCounts();
  
public:
  // constructors
  NavierStokesSolver(){ };
  NavierStokesSolver(CartesianMesh *cartesianMesh, 
                     FlowDescription<dim> *flowDescription, 
                     SimulationParameters *simulationParameters);
  // destructor
  ~NavierStokesSolver(){ };

  // initialize systems
  virtual PetscErrorCode initialize();
  // clean up at end of simulation
  virtual PetscErrorCode finalize();

  // evaluate if the simulation is completed
  PetscBool finished();

  // write numerical solution into respective files
  virtual PetscErrorCode writeData();
  // read numerical solution from files
  PetscErrorCode readData();

}; // NavierStokesSolver

#endif