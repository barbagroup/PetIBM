/*! Implementation of the methods of the class `LiEtAlSolver`.
 * \file LiEtAlSolver.cpp
 */


#include "LiEtAlSolver.h"

#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <sys/stat.h>

#include <petscdmcomposite.h>


/*!
 * \brief Constructor.
 *
 * The constructor calls the `NavierStokesSolver` constructor and initializes
 * additional pointers related to the Lagrangian points.
 */
template <PetscInt dim>
LiEtAlSolver<dim>::LiEtAlSolver(CartesianMesh *cartesianMesh, 
                                FlowDescription<dim> *flowDescription, 
                                SimulationParameters *simulationParameters)
                    : NavierStokesSolver<dim>::NavierStokesSolver(cartesianMesh, 
                                                                  flowDescription, 
                                                                  simulationParameters)
{
  fTilde = PETSC_NULL;
  rhsf = PETSC_NULL;
  temp = PETSC_NULL;
  temp2 = PETSC_NULL;
  E = PETSC_NULL;
  ET = PETSC_NULL;
  EBNET = PETSC_NULL;
  PetscLogStageRegister("RHSForce", &stageRHSForceSystem);
  PetscLogStageRegister("solveForce", &stageSolveForceSystem);
} // LiEtAl


/*!
 * \brief Initializes the solver.
 */
template <PetscInt dim>
PetscErrorCode LiEtAlSolver<dim>::initialize()
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  ierr = PetscLogStagePush(NavierStokesSolver<dim>::stageInitialize); CHKERRQ(ierr);

  ierr = initializeBodies(); CHKERRQ(ierr);
  ierr = createDMs(); CHKERRQ(ierr);
  ierr = createVecs(); CHKERRQ(ierr);
  ierr = createGlobalMappingBodies(); CHKERRQ(ierr);

  ierr = NavierStokesSolver<dim>::initializeCommon(); CHKERRQ(ierr);

  ierr = generateET(); CHKERRQ(ierr);
  ierr = MatTranspose(ET, MAT_INITIAL_MATRIX, &E); CHKERRQ(ierr);
  ierr = MatDiagonalScale(ET, NavierStokesSolver<dim>::BN, NULL); CHKERRQ(ierr);
  ierr = MatMatMult(E, ET, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &EBNET); CHKERRQ(ierr);
  ierr = PetscObjectViewFromOptions((PetscObject) EBNET, NULL, "-EBNET_mat_view"); CHKERRQ(ierr);

  ierr = createForceSolver(); CHKERRQ(ierr);

  ierr = PetscLogStagePop(); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // initialize


/**
 * \brief Advance in time. Calculates the variables at the next time-step.
 */
template <PetscInt dim>
PetscErrorCode LiEtAlSolver<dim>::stepTime()
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  if (dim == 2)
  {
    ierr = DMCompositeScatter(NavierStokesSolver<dim>::qPack,
                              NavierStokesSolver<dim>::q,
                              NavierStokesSolver<dim>::qxLocal,
                              NavierStokesSolver<dim>::qyLocal); CHKERRQ(ierr);
  }
  else
  {
    ierr = DMCompositeScatter(NavierStokesSolver<dim>::qPack,
                              NavierStokesSolver<dim>::q,
                              NavierStokesSolver<dim>::qxLocal,
                              NavierStokesSolver<dim>::qyLocal,
                              NavierStokesSolver<dim>::qzLocal); CHKERRQ(ierr);
  }

  // solve system for intermediate velocity
  ierr = NavierStokesSolver<dim>::assembleRHSVelocity(); CHKERRQ(ierr);
  ierr = NavierStokesSolver<dim>::solveIntermediateVelocity(); CHKERRQ(ierr);

  PetscInt maxIters = 5;
  for (PetscInt iter=0; iter<maxIters; iter++)
  {
    if (iter > 0)
    {
      ierr = VecCopy(NavierStokesSolver<dim>::q, NavierStokesSolver<dim>::qStar); CHKERRQ(ierr);
    }

    // solve system for Lagrangian forces
    ierr = assembleRHSForce(); CHKERRQ(ierr);
    ierr = solveForceSystem(); CHKERRQ(ierr);
    ierr = updateFlux();

    // solve Poisson system for the pressure
    ierr = NavierStokesSolver<dim>::assembleRHSPoisson(); CHKERRQ(ierr);
    ierr = NavierStokesSolver<dim>::solvePoissonSystem(); CHKERRQ(ierr);

    // project intermediate velocity field to satisfy divergence-free condition
    ierr = NavierStokesSolver<dim>::projectionStep(); CHKERRQ(ierr);
  }

  NavierStokesSolver<dim>::timeStep++;

  PetscFunctionReturn(0);
} // stepTime


/**
 * \brief Assembles the right-hand side of the system for the Lagrangian forces.
 */
template <PetscInt dim>
PetscErrorCode LiEtAlSolver<dim>::assembleRHSForce()
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  ierr = PetscLogStagePush(stageRHSForceSystem); CHKERRQ(ierr);

  ierr = MatMult(E, NavierStokesSolver<dim>::qStar, rhsf);

  ierr = PetscObjectViewFromOptions((PetscObject) rhsf, NULL, "-rhsf_vec_view"); CHKERRQ(ierr);

  ierr = PetscLogStagePop(); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // assembleRHSForce


/**
 * \brief Solves the system for the Lagrangian forces.
 */
template <PetscInt dim>
PetscErrorCode LiEtAlSolver<dim>::solveForceSystem()
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  ierr = PetscLogStagePush(stageSolveForceSystem); CHKERRQ(ierr);

  ierr = forces->solve(fTilde, rhsf); CHKERRQ(ierr);
  
  ierr = PetscLogStagePop(); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // solveForceSystem


/**
 * \brief Updates the intermediate fluxes to satisfy the no-slip condition.
 */
template <PetscInt dim>
PetscErrorCode LiEtAlSolver<dim>::updateFlux()
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  ierr = PetscLogStagePush(NavierStokesSolver<dim>::stageProjectionStep); CHKERRQ(ierr);

  ierr = MatMult(ET, fTilde, temp); CHKERRQ(ierr);
  ierr = VecPointwiseMult(temp2, NavierStokesSolver<dim>::BN, temp); CHKERRQ(ierr);
  ierr = VecCopy(NavierStokesSolver<dim>::qStar, qStar2); CHKERRQ(ierr);
  ierr = VecWAXPY(NavierStokesSolver<dim>::qStar, -1.0, temp2, qStar2); CHKERRQ(ierr);

  ierr = PetscLogStagePop(); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // updateFlux


/*!
 * \brief Gets the total number of Lagrangian points.
 *
 * \param n The integer to return.
 */
template <PetscInt dim>
PetscErrorCode LiEtAlSolver<dim>::getNumLagPoints(PetscInt &n)
{
  PetscFunctionBeginUser;

  n = 0;
  for (auto &body : bodies)
  {
    n += body.numPoints;
  }
  
  PetscFunctionReturn(0);
} // getNumLagPoints


/*!
 * \brief Gets the number of Lagrangian points on each process.
 *
 * \param numOnProcess The vector of integers to return.
 */
template <PetscInt dim>
PetscErrorCode LiEtAlSolver<dim>::getNumLagPointsOnProcess(std::vector<PetscInt> &numOnProcess)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  PetscMPIInt rank, size;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);

  numOnProcess.resize(size);
  numOnProcess[rank] = 0;
  for (auto &body : bodies)
  {
    numOnProcess[rank] += body.numPointsOnProcess[rank];
  }

  ierr = MPI_Allgather(MPI_IN_PLACE, 1, MPIU_INT,
                       &numOnProcess[0], 1, MPIU_INT,
                       PETSC_COMM_WORLD); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // getNumLagPointsOnProcess


/*!
 * \brief Registers the Lagrangian points on local process for each body.
 *
 * A local process is defined by a box (xmin, xmax, ymin, ymax) whose dimensions
 * are set using the layout from the DMDA object for the pressure.
 */
template <PetscInt dim>
PetscErrorCode LiEtAlSolver<dim>::registerLagPointsOnProcess()
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  DMDALocalInfo info;
  ierr = DMDAGetLocalInfo(NavierStokesSolver<dim>::pda, &info); CHKERRQ(ierr);
  PetscReal box[2*dim];
  box[0] = NavierStokesSolver<dim>::mesh->x[info.xs];
  box[1] = NavierStokesSolver<dim>::mesh->x[info.xs + info.xm];
  box[2] = NavierStokesSolver<dim>::mesh->y[info.ys];
  box[3] = NavierStokesSolver<dim>::mesh->y[info.ys + info.ym];
  if (dim == 3)
  {
    box[4] = NavierStokesSolver<dim>::mesh->z[info.zs];
    box[5] = NavierStokesSolver<dim>::mesh->z[info.zs + info.zm];
  }

  for (auto &body : bodies)
  {
    ierr = body.registerPointsOnProcess(box); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
} // registerLagPointsOnProcess


/**
 * \brief Destroys PETSc objects.
 */
template <PetscInt dim>
PetscErrorCode LiEtAlSolver<dim>::finalize()
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  ierr = NavierStokesSolver<dim>::finalize();
  if (bda != PETSC_NULL) {ierr = DMDestroy(&bda); CHKERRQ(ierr);}
  if (fTilde != PETSC_NULL) {ierr = VecDestroy(&fTilde); CHKERRQ(ierr);}
  if (rhsf != PETSC_NULL) {ierr = VecDestroy(&rhsf); CHKERRQ(ierr);}
  if (qStar2 != PETSC_NULL) {ierr = VecDestroy(&qStar2); CHKERRQ(ierr);}
  if (temp != PETSC_NULL) {ierr = VecDestroy(&temp); CHKERRQ(ierr);}
  if (temp2 != PETSC_NULL) {ierr = VecDestroy(&temp2); CHKERRQ(ierr);}
  if (E != PETSC_NULL) {ierr = MatDestroy(&E); CHKERRQ(ierr);}
  if (ET != PETSC_NULL) {ierr = MatDestroy(&ET); CHKERRQ(ierr);}
  if (EBNET != PETSC_NULL) {ierr = MatDestroy(&EBNET); CHKERRQ(ierr);}
  
  delete forces;

  PetscFunctionReturn(0);
}  // finalize


#include "inline/initializeBodies.inl"
#include "inline/createDMs.inl"
#include "inline/createGlobalMappingBodies.inl"
#include "inline/createVecs.inl"
#include "inline/generateET.inl"
#include "inline/createForceSolver.inl"
#include "inline/calculateForces.inl"
#include "inline/io.inl"


// dimensions specialization
template class LiEtAlSolver<2>;
template class LiEtAlSolver<3>;
