/***************************************************************************//**
 * \file NavierStokesSolver.cpp
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the methods of the class `NavierStokesSolver`.
 */


#include "NavierStokesSolver.h"

#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <sys/stat.h>

#include <petscdmcomposite.h>


/**
 * \brief Constructor: Stores simulation parameters and initializes pointers.
 *
 * \param cartesianMesh Contains info about the Cartesian grid.
 * \param flowDescription Contains info about the fluid and the flow.
 * \param simulationParameters Contains info about the parameters of the simulation.
 */
template <PetscInt dim>
NavierStokesSolver<dim>::NavierStokesSolver(CartesianMesh *cartesianMesh, 
                                            FlowDescription<dim> *flowDescription, 
                                            SimulationParameters *simulationParameters) 
{
  // simulation info
  mesh = cartesianMesh;
  flow = flowDescription;
  parameters = simulationParameters;
  timeStep = parameters->startStep;
  // DM objects
  lambdaPack = PETSC_NULL;
  pda = PETSC_NULL;
  qPack = PETSC_NULL;
  uda = PETSC_NULL;
  vda = PETSC_NULL;
  wda = PETSC_NULL;
  // global vectors
  q = PETSC_NULL;
  qStar = PETSC_NULL;
  H = PETSC_NULL;
  rn = PETSC_NULL;
  bc1 = PETSC_NULL;
  rhs1 = PETSC_NULL;
  lambda = PETSC_NULL;
  r2 = PETSC_NULL;
  rhs2 = PETSC_NULL;
  temp = PETSC_NULL;
  // local vectors
  qxLocal = PETSC_NULL;
  qyLocal = PETSC_NULL;
  qzLocal = PETSC_NULL;
  // mappings local to global
  pMapping = PETSC_NULL;
  uMapping = PETSC_NULL;
  vMapping = PETSC_NULL;
  wMapping = PETSC_NULL;
  // Mats
  A = PETSC_NULL;
  QT = PETSC_NULL;
  BNQ = PETSC_NULL;
  QTBNQ = PETSC_NULL;
  // diagonal matrices
  BN = PETSC_NULL;
  RInv = PETSC_NULL;
  MHat = PETSC_NULL;
  // KSPs
  ksp1 = PETSC_NULL;
  ksp2 = PETSC_NULL;
  // PetscLogStages
  PetscLogStageRegister("initialize", &stageInitialize);
  PetscLogStageRegister("RHSVelocity", &stageRHSVelocitySystem);
  PetscLogStageRegister("solveVelocity", &stageSolveVelocitySystem);
  PetscLogStageRegister("RHSPoisson", &stageRHSPoissonSystem);
  PetscLogStageRegister("solvePoisson", &stageSolvePoissonSystem);
  PetscLogStageRegister("projectionStep", &stageProjectionStep);
} // NavierStokesSolver


/**
 * \brief Initializes the solver.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::initialize()
{
  PetscErrorCode ierr;

  ierr = PetscLogStagePush(stageInitialize); CHKERRQ(ierr);
  
  ierr = printInfo(); CHKERRQ(ierr);
  ierr = createDMs(); CHKERRQ(ierr);
  ierr = initializeCommon(); CHKERRQ(ierr);
  ierr = writeGrid(); CHKERRQ(ierr);
  
  ierr = PetscLogStagePop(); CHKERRQ(ierr);

  return 0;
} // initialize


/**
 * \brief Initializes data common to `NavierStokesSolver` and its derived classes.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::initializeCommon()
{
  PetscErrorCode ierr;

  ierr = createVecs(); CHKERRQ(ierr);
  
  ierr = initializeFluxes(); CHKERRQ(ierr);
  ierr = initializeLambda(); CHKERRQ(ierr);
  ierr = updateBoundaryGhosts(); CHKERRQ(ierr);

  ierr = createLocalToGlobalMappingsFluxes(); CHKERRQ(ierr);
  ierr = createLocalToGlobalMappingsLambda(); CHKERRQ(ierr);

  ierr = generateDiagonalMatrices(); CHKERRQ(ierr);
  ierr = generateA(); CHKERRQ(ierr);
  ierr = generateBNQ(); CHKERRQ(ierr);
  ierr = generateQTBNQ(); CHKERRQ(ierr);
  ierr = createKSPs(); CHKERRQ(ierr);
  ierr = setNullSpace(); CHKERRQ(ierr);

  return 0;
} // initializeCommon


/**
 * \brief Adavance in time. Calculates the variables at the next time-step.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::stepTime()
{
  PetscErrorCode ierr;

  // solve system for intermediate velocity
  ierr = assembleRHSVelocity(); CHKERRQ(ierr);
  ierr = solveIntermediateVelocity(); CHKERRQ(ierr);

  // solve Poisson system for Lagrange multipliers
  // Navier-Stokes: Lagrange multipliers = pressure
  // Taira-Colonius: Lagrange multipliers = pressure + body forces
  ierr = assembleRHSPoisson(); CHKERRQ(ierr);
  ierr = solvePoissonSystem(); CHKERRQ(ierr);

  // project intermediate velocity field to satisfy divergence-free condition
  // and no-slip condition at immersed boundary (when Taira-Colonius method used)
  ierr = projectionStep(); CHKERRQ(ierr);
  
  // code-development helpers: output vectors and matrices
  ierr = helpers(); CHKERRQ(ierr);

  timeStep++;

  return 0;
} // stepTime


/**
 * \brief Assembles the right hand-side of the velocity system.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::assembleRHSVelocity()
{
  PetscErrorCode ierr;

  ierr = PetscLogStagePush(stageRHSVelocitySystem); CHKERRQ(ierr);

  ierr = calculateExplicitTerms(); CHKERRQ(ierr);
  ierr = updateBoundaryGhosts(); CHKERRQ(ierr);
  ierr = generateBC1(); CHKERRQ(ierr);
  ierr = VecWAXPY(rhs1, 1.0, rn, bc1); CHKERRQ(ierr);
  ierr = VecPointwiseMult(rhs1, MHat, rhs1); CHKERRQ(ierr);

  ierr = PetscLogStagePop(); CHKERRQ(ierr);

  return 0;
} // assembleRHSVelocity


/**
 * \brief Solves the system for the intermediate fluxes.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::solveIntermediateVelocity()
{
  PetscErrorCode ierr;

  ierr = PetscLogStagePush(stageSolveVelocitySystem); CHKERRQ(ierr);

  ierr = KSPSolve(ksp1, rhs1, qStar); CHKERRQ(ierr);
  
  ierr = PetscLogStagePop(); CHKERRQ(ierr);

  KSPConvergedReason reason;
  ierr = KSPGetConvergedReason(ksp1, &reason); CHKERRQ(ierr);
  if (reason < 0)
  {
    ierr = PetscPrintf(PETSC_COMM_WORLD, "\n[time-step %d]", timeStep); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,
                       "\nERROR: velocity solver diverged due to reason: %d\n", 
                       reason); CHKERRQ(ierr);
    exit(0);
  }

  return 0;
} // solveIntermediateVelocity


/**
 * \brief Assembles the right hand-side of the Poisson system.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::assembleRHSPoisson()
{
  PetscErrorCode ierr;

  ierr = PetscLogStagePush(stageRHSPoissonSystem); CHKERRQ(ierr);

  ierr = generateR2(); CHKERRQ(ierr);
  ierr = VecScale(r2, -1.0); CHKERRQ(ierr);
  ierr = MatMultAdd(QT, qStar, r2, rhs2); CHKERRQ(ierr);
  
  ierr = PetscLogStagePop(); CHKERRQ(ierr);  

  return 0;
} // assembleRHSPoisson


/**
 * \brief Solves Poisson system for the pressure-forces.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::solvePoissonSystem()
{
  PetscErrorCode ierr;
  
  ierr = PetscLogStagePush(stageSolvePoissonSystem); CHKERRQ(ierr);

  ierr = KSPSolve(ksp2, rhs2, lambda); CHKERRQ(ierr);
  
  ierr = PetscLogStagePop(); CHKERRQ(ierr);
  
  KSPConvergedReason reason;
  ierr = KSPGetConvergedReason(ksp2, &reason); CHKERRQ(ierr);
  if (reason < 0)
  {
    ierr = PetscPrintf(PETSC_COMM_WORLD, "\n[time-step %d]", timeStep); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,
                       "\nERROR: Poisson solver diverged due to reason: %d\n", 
                       reason); CHKERRQ(ierr);
    exit(0);
  }

  return 0;
} // solvePoissonSystem


/**
 * \brief Projects the fluxes onto the space satisfying the constraints.
 *
 * In the case of a pure Navier-Stokes solver, the intermediate velocity is 
 * projected on the divergence-free space.
 * In the case where the immersed boundary projection method is used, the 
 * velocity field is also projected onto the space where the no-slip condition 
 * is satisfied.
 *
 * \f[ q = q^* - B^N Q \lambda \f]
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::projectionStep()
{
  PetscErrorCode ierr;

  ierr = PetscLogStagePush(stageProjectionStep); CHKERRQ(ierr);
  
  ierr = MatMult(BNQ, lambda, temp); CHKERRQ(ierr);
  ierr = VecWAXPY(q, -1.0, temp, qStar); CHKERRQ(ierr);

  ierr = PetscLogStagePop(); CHKERRQ(ierr);

  return 0;
} // projectionStep


/**
 * \brief Computes the matrix \f$ Q^T B^N Q \f$.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::generateQTBNQ()
{
  PetscErrorCode ierr;

  PetscLogEvent GENERATE_QTBNQ;
  ierr = PetscLogEventRegister("generateQTBNQ", 0, &GENERATE_QTBNQ); CHKERRQ(ierr);
  ierr = PetscLogEventBegin(GENERATE_QTBNQ, 0, 0, 0, 0); CHKERRQ(ierr);

  ierr = MatMatMult(QT, BNQ, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &QTBNQ); CHKERRQ(ierr);
  
  ierr = PetscLogEventEnd(GENERATE_QTBNQ, 0, 0, 0, 0); CHKERRQ(ierr);

  return 0;
} // generateQTBNQ


/**
 * \brief Count the numbers of non-zeros in the diagonal 
 *        and off-diagonal portions of the parallel matrices.
 *
 * \param cols Array of column indices where non-zeros are present in a particular
 *             row of a matrix
 * \param numCols Number of column indices in the given array
 * \param rowStart Starting index of the portion of the result vector 
 *                 that resides on the current process
 * \param rowEnd Ending index, which is 1 greater than the index of the last row
 *               of the portion of the result vector that resides on the current
 *               process
 * \param d_nnz Number of non-zeros on the diagonal portion of the matrix
 * \param o_nnz Number of non-zeros off-diagonal of the portion of the matrix
 *
 * `d_nnz` and `o_nnz` are passed by reference, and are outputs of the function.
 *
 */
template <PetscInt dim>
void NavierStokesSolver<dim>::countNumNonZeros(PetscInt *cols, size_t numCols, 
                                               PetscInt rowStart, PetscInt rowEnd, 
                                               PetscInt &d_nnz, PetscInt &o_nnz)
{
  d_nnz = 0;
  o_nnz = 0;
  for(size_t i=0; i<numCols; i++)
  {
    (cols[i]>=rowStart && cols[i]<rowEnd)? d_nnz++ : o_nnz++;
  }
} // countNumNonZeros


/**
 * \brief Is the simulation completed?
 */
template <PetscInt dim>
PetscBool NavierStokesSolver<dim>::finished()
{
  return (timeStep >= parameters->startStep+parameters->nt)? PETSC_TRUE : PETSC_FALSE;
} // finished


/**
 * \brief Code-development helpers: outputs vectors and matrices to files.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::helpers()
{
  PetscErrorCode ierr;

  if (timeStep == parameters->startStep+1)
  {
    PetscBool outputToFiles = PETSC_FALSE;

    ierr = PetscOptionsBegin(PETSC_COMM_WORLD, nullptr, nullptr, nullptr); CHKERRQ(ierr);
    ierr = PetscOptionsBool("-outputs", 
                            "Outputs vectors and matrices to check implementation", "",
                            outputToFiles, &outputToFiles, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);
    if (outputToFiles)
    {
      ierr = helperOutputVectors(); CHKERRQ(ierr);
      ierr = helperOutputMatrices(); CHKERRQ(ierr);
      ierr = finalize();
      ierr = PetscFinalize(); CHKERRQ(ierr);
      exit(0);
    }
  }

  return 0;
} // helpers


/**
 * \brief Deallocate memory to avoid memory leaks.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::finalize()
{
  PetscErrorCode ierr;
  
  // DM objects
  if (pda != PETSC_NULL)       {ierr = DMDestroy(&pda); CHKERRQ(ierr);}
  if (uda != PETSC_NULL)       {ierr = DMDestroy(&uda); CHKERRQ(ierr);}
  if (vda != PETSC_NULL)       {ierr = DMDestroy(&vda); CHKERRQ(ierr);}
  if (wda != PETSC_NULL)       {ierr = DMDestroy(&wda); CHKERRQ(ierr);}
  if (qPack != PETSC_NULL)     {ierr = DMDestroy(&qPack); CHKERRQ(ierr);}
  if (lambdaPack != PETSC_NULL){ierr = DMDestroy(&lambdaPack); CHKERRQ(ierr);}
  // global solution vectors
  if (q != PETSC_NULL)    {ierr = VecDestroy(&q); CHKERRQ(ierr);}
  if (qStar != PETSC_NULL){ierr = VecDestroy(&qStar); CHKERRQ(ierr);}
  // local fluxes vectors
  if (qxLocal != PETSC_NULL){ierr = VecDestroy(&qxLocal); CHKERRQ(ierr);}
  if (qyLocal != PETSC_NULL){ierr = VecDestroy(&qyLocal); CHKERRQ(ierr);}
  if (qzLocal != PETSC_NULL){ierr = VecDestroy(&qzLocal); CHKERRQ(ierr);}
  // global vectors
  if (H != PETSC_NULL)     {ierr = VecDestroy(&H); CHKERRQ(ierr);}
  if (rn != PETSC_NULL)    {ierr = VecDestroy(&rn); CHKERRQ(ierr);}
  if (bc1 != PETSC_NULL)   {ierr = VecDestroy(&bc1); CHKERRQ(ierr);}
  if (rhs1 != PETSC_NULL)  {ierr = VecDestroy(&rhs1); CHKERRQ(ierr);}
  if (temp != PETSC_NULL)  {ierr = VecDestroy(&temp); CHKERRQ(ierr);}
  if (lambda != PETSC_NULL){ierr = VecDestroy(&lambda); CHKERRQ(ierr);}
  if (r2 != PETSC_NULL)    {ierr = VecDestroy(&r2); CHKERRQ(ierr);}
  if (rhs2 != PETSC_NULL)  {ierr = VecDestroy(&rhs2); CHKERRQ(ierr);}
  // mappings local to global indices
  if (uMapping != PETSC_NULL){ierr = VecDestroy(&uMapping); CHKERRQ(ierr);}
  if (vMapping != PETSC_NULL){ierr = VecDestroy(&vMapping); CHKERRQ(ierr);}
  if (wMapping != PETSC_NULL){ierr = VecDestroy(&wMapping); CHKERRQ(ierr);}
  if (pMapping != PETSC_NULL){ierr = VecDestroy(&pMapping); CHKERRQ(ierr);}
  // diagonal matrices
  if (MHat != PETSC_NULL){ierr = VecDestroy(&MHat); CHKERRQ(ierr);}
  if (RInv != PETSC_NULL){ierr = VecDestroy(&RInv); CHKERRQ(ierr);}
  if (BN != PETSC_NULL)  {ierr = VecDestroy(&BN); CHKERRQ(ierr);}
  // matrices
  if (A != PETSC_NULL)    {ierr = MatDestroy(&A); CHKERRQ(ierr);}
  if (QT != PETSC_NULL)   {ierr = MatDestroy(&QT); CHKERRQ(ierr);}
  if (BNQ != PETSC_NULL)  {ierr = MatDestroy(&BNQ); CHKERRQ(ierr);}
  if (QTBNQ != PETSC_NULL){ierr = MatDestroy(&QTBNQ); CHKERRQ(ierr);}
  // KSPs
  if (ksp1 != PETSC_NULL){ierr = KSPDestroy(&ksp1); CHKERRQ(ierr);}
  if (ksp2 != PETSC_NULL){ierr = KSPDestroy(&ksp2); CHKERRQ(ierr);}

  // print performance summary to file
  PetscViewer viewer;
  std::string filePath = parameters->directory + "/performanceSummary.txt";
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, filePath.c_str(), &viewer); CHKERRQ(ierr);
  ierr = PetscLogView(viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  return 0;
} // finalize


#include "inline/createDMs.inl"
#include "inline/createVecs.inl"
#include "inline/createKSPs.inl"
#include "inline/setNullSpace.inl"
#include "inline/createLocalToGlobalMappingsFluxes.inl"
#include "inline/createLocalToGlobalMappingsLambda.inl"
#include "inline/addInitialPerturbation.inl"
#include "inline/initializeFluxes.inl"
#include "inline/initializeLambda.inl"
#include "inline/updateBoundaryGhosts.inl"
#include "inline/calculateExplicitTerms.inl"
#include "inline/generateDiagonalMatrices.inl"
#include "inline/generateA.inl"
#include "inline/generateBC1.inl"
#include "inline/generateBNQ.inl"
#include "inline/generateR2.inl"
#include "inline/io.inl"


// template class specialization
template class NavierStokesSolver<2>;
template class NavierStokesSolver<3>;
