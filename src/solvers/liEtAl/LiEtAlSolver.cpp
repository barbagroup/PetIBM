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
  dfTilde = PETSC_NULL;
  rhsf = PETSC_NULL;
  tmp = PETSC_NULL;
  dlambda = PETSC_NULL;
  rhs1_n = PETSC_NULL;
  E = PETSC_NULL;
  ET = PETSC_NULL;
  EBNET = PETSC_NULL;
  G = PETSC_NULL;
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

  ierr = NavierStokesSolver<dim>::generateGradient(&G); CHKERRQ(ierr);
  ierr = generateET(); CHKERRQ(ierr);
  ierr = MatTranspose(ET, MAT_INITIAL_MATRIX, &E); CHKERRQ(ierr);
  Mat BNET;
  ierr = MatDuplicate(ET, MAT_COPY_VALUES, &BNET);
  ierr = MatDiagonalScale(BNET, NavierStokesSolver<dim>::BN, NULL); CHKERRQ(ierr);
  ierr = MatMatMult(E, BNET, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &EBNET); CHKERRQ(ierr);
  ierr = MatDestroy(&BNET); CHKERRQ(ierr);
  ierr = PetscObjectViewFromOptions((PetscObject) EBNET, NULL, "-EBNET_mat_view"); CHKERRQ(ierr);

  ierr = createForceSolver(); CHKERRQ(ierr);

  // info about the algorithm to use and the inner-iterations
  algorithm = NavierStokesSolver<dim>::parameters->decoupling_algorithm;
  forceEstimator = NavierStokesSolver<dim>::parameters->decoupling_forceEstimator;
  atol = NavierStokesSolver<dim>::parameters->decoupling_atol;
  rtol = NavierStokesSolver<dim>::parameters->decoupling_rtol;
  maxIters = NavierStokesSolver<dim>::parameters->decoupling_maxIters;
  printStats = NavierStokesSolver<dim>::parameters->decoupling_printStats;

  ierr = PetscLogStagePop(); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // initialize


/*!
 * \brief Advance in time. Calculates the variables at the next time-step.
 */
template <PetscInt dim>
PetscErrorCode LiEtAlSolver<dim>::stepTime()
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  ierr = scatterGlobalToLocal(); CHKERRQ(ierr);
  
  // estimation of the momentum forcing at beginning of time step
  // forceEstimator = 2: use forcing from previous time step
  if (forceEstimator == 1)
  {
    ierr = VecSet(fTilde, 0.0); CHKERRQ(ierr);
  }
  else if (forceEstimator == 3)
  {
    ierr = assembleRHSForce(NavierStokesSolver<dim>::q); CHKERRQ(ierr);
    ierr = solveForceSystem(fTilde); CHKERRQ(ierr);
  }

  ierr = NavierStokesSolver<dim>::assembleRHSVelocity(); CHKERRQ(ierr);
  if (maxIters > 1)
  {
    // hold the convective terms and diffusive terms
    ierr = VecCopy(NavierStokesSolver<dim>::rhs1, rhs1_n); CHKERRQ(ierr);
  }

  PetscInt iter = 0;
  PetscReal norm = 1.0, norm_init = 1.0, ratio;
  while (norm > std::max(atol, rtol * norm_init) && iter < maxIters)
  {
    ierr = updateRHSVelocity(); CHKERRQ(ierr);
    ierr = NavierStokesSolver<dim>::solveIntermediateVelocity(); CHKERRQ(ierr);

    if (algorithm == 1)
    {
      ierr = assembleRHSForce(NavierStokesSolver<dim>::qStar); CHKERRQ(ierr);
      ierr = solveForceSystem(dfTilde); CHKERRQ(ierr);
      ierr = updateFlux(dfTilde);

      ierr = NavierStokesSolver<dim>::assembleRHSPoisson(); CHKERRQ(ierr);
      ierr = solvePoissonSystem(dlambda); CHKERRQ(ierr);
      ierr = projectionStep(dlambda); CHKERRQ(ierr);
    }
    else if (algorithm == 3)
    {
      ierr = NavierStokesSolver<dim>::assembleRHSPoisson(); CHKERRQ(ierr);
      ierr = solvePoissonSystem(dlambda); CHKERRQ(ierr);
      ierr = projectionStep(dlambda); CHKERRQ(ierr);

      ierr = assembleRHSForce(NavierStokesSolver<dim>::qStar); CHKERRQ(ierr);
      ierr = solveForceSystem(dfTilde); CHKERRQ(ierr);
      ierr = updateFlux(dfTilde);
    }
    else
    {
      SETERRQ1(PETSC_COMM_WORLD, 63,
             "Algorithm id %d not recognized (values accepted: 1 and 3)",
             algorithm);
    }

    // update Lagrangian forces and pressure field
    ierr = VecAXPY(fTilde, 1.0, dfTilde); CHKERRQ(ierr);
    ierr = VecAXPY(NavierStokesSolver<dim>::lambda, 1.0, dlambda); CHKERRQ(ierr);

    if (maxIters > 1)
    {
      ierr = VecNorm(dfTilde, NORM_2, &norm); CHKERRQ(ierr);
      ierr = VecNorm(fTilde, NORM_2, &norm_init); CHKERRQ(ierr);
      ratio = norm / norm_init;
      if (printStats)
      {
        ierr = PetscPrintf(PETSC_COMM_WORLD,
                           "[time-step %d][iter %d] L2(df_k)=%f\tL2(df_k)/L2(f_k)=%f\n",
                           NavierStokesSolver<dim>::timeStep, iter, norm, ratio); CHKERRQ(ierr);
      }
    }
    iter++;
  }

  NavierStokesSolver<dim>::timeStep++;

  PetscFunctionReturn(0);
} // stepTime


/*!
 * \brief Scatter global flux vector into local flux vectors.
 */
template <PetscInt dim>
PetscErrorCode LiEtAlSolver<dim>::scatterGlobalToLocal()
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

  PetscFunctionReturn(0);
} // scatterGlobalToLocal


/*!
 * \brief Assembles the right hand-side of the velocity system.
 */
template <PetscInt dim>
PetscErrorCode LiEtAlSolver<dim>::updateRHSVelocity()
{
  PetscErrorCode ierr;

  ierr = PetscLogStagePush(NavierStokesSolver<dim>::stageRHSVelocitySystem); CHKERRQ(ierr);

  if (maxIters > 1)
  {
    ierr = VecCopy(rhs1_n, NavierStokesSolver<dim>::rhs1); CHKERRQ(ierr);
  }
  // add estimation of the momentum forcing
  ierr = MatMult(ET, fTilde, tmp); CHKERRQ(ierr);
  ierr = VecAXPY(NavierStokesSolver<dim>::rhs1, -1.0, tmp); CHKERRQ(ierr);
  // add estimation of the pressure gradient
  ierr = MatMult(G, NavierStokesSolver<dim>::lambda, tmp); CHKERRQ(ierr);
  ierr = VecAXPY(NavierStokesSolver<dim>::rhs1, -1.0, tmp); CHKERRQ(ierr);
  ierr = PetscObjectViewFromOptions((PetscObject) NavierStokesSolver<dim>::rhs1, NULL, "-rhs1_vec_view"); CHKERRQ(ierr);

  ierr = PetscLogStagePop(); CHKERRQ(ierr);

  return 0;
} // updateRHSVelocity


/*!
 * \brief Assembles the right-hand side of the system for the Lagrangian forces.
 *
 * \param q Flux vector to use in the RHS.
 */
template <PetscInt dim>
PetscErrorCode LiEtAlSolver<dim>::assembleRHSForce(Vec q)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  ierr = PetscLogStagePush(stageRHSForceSystem); CHKERRQ(ierr);

  ierr = MatMult(E, q, rhsf);

  ierr = PetscObjectViewFromOptions((PetscObject) rhsf, NULL, "-rhsf_vec_view"); CHKERRQ(ierr);

  ierr = PetscLogStagePop(); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // assembleRHSForce


/*!
 * \brief Solves the Lagrangian system.
 *
 * \param f The solution vector.
 */
template <PetscInt dim>
PetscErrorCode LiEtAlSolver<dim>::solveForceSystem(Vec &f)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  ierr = PetscLogStagePush(stageSolveForceSystem); CHKERRQ(ierr);

  ierr = forces->solve(f, rhsf); CHKERRQ(ierr);
  
  ierr = PetscLogStagePop(); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // solveForceSystem


/*!
 * \brief Updates the intermediate fluxes to satisfy the no-slip condition.
 *
 * \param f The force vector or the delta-force vector.
 */
template <PetscInt dim>
PetscErrorCode LiEtAlSolver<dim>::updateFlux(Vec f)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  ierr = PetscLogStagePush(NavierStokesSolver<dim>::stageProjectionStep); CHKERRQ(ierr);

  ierr = MatMult(ET, f, tmp); CHKERRQ(ierr);
  ierr = VecPointwiseMult(tmp, NavierStokesSolver<dim>::BN, tmp); CHKERRQ(ierr);
  if (algorithm == 1)
  {
    ierr = VecAXPY(NavierStokesSolver<dim>::qStar, -1.0, tmp); CHKERRQ(ierr);  
  }
  else if (algorithm == 3)
  {
    ierr = VecWAXPY(NavierStokesSolver<dim>::q, -1.0, tmp, NavierStokesSolver<dim>::qStar); CHKERRQ(ierr);
  }
  else
  {
    SETERRQ1(PETSC_COMM_WORLD, 63,
             "Algorithm id %d not recognized (values accepted: 1 and 3)",
             algorithm);
  }

  ierr = PetscLogStagePop(); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // updateFlux


/*!
 * \brief Solves Poisson system.
 *
 * \param p The solution vector.
 */
template <PetscInt dim>
PetscErrorCode LiEtAlSolver<dim>::solvePoissonSystem(Vec &p)
{
  PetscErrorCode ierr;
  
  ierr = PetscLogStagePush(NavierStokesSolver<dim>::stageSolvePoissonSystem); CHKERRQ(ierr);

  ierr = NavierStokesSolver<dim>::poisson->solve(p, NavierStokesSolver<dim>::rhs2); CHKERRQ(ierr);
  
  ierr = PetscLogStagePop(); CHKERRQ(ierr);

  return 0;
} // solvePoissonSystem


/*!
 * \brief Project fluxes on divergence-free space.
 *
 * \param p The pressure vector and the delta-pressure vector.
 */
template <PetscInt dim>
PetscErrorCode LiEtAlSolver<dim>::projectionStep(Vec p)
{
  PetscErrorCode ierr;

  ierr = PetscLogStagePush(NavierStokesSolver<dim>::stageProjectionStep); CHKERRQ(ierr);
  
  ierr = MatMult(NavierStokesSolver<dim>::BNQ, p, NavierStokesSolver<dim>::temp); CHKERRQ(ierr);
  if (algorithm == 1)
  {
    ierr = VecWAXPY(NavierStokesSolver<dim>::q, -1.0, NavierStokesSolver<dim>::temp, NavierStokesSolver<dim>::qStar); CHKERRQ(ierr);  
  }
  else if (algorithm == 3)
  {
    ierr = VecAXPY(NavierStokesSolver<dim>::qStar, -1.0, NavierStokesSolver<dim>::temp); CHKERRQ(ierr);
  }
  else
  {
    SETERRQ1(PETSC_COMM_WORLD, 63,
             "Algorithm id %d not recognized (values accepted: 1 and 3)",
             algorithm);
  }

  ierr = PetscLogStagePop(); CHKERRQ(ierr);

  return 0;
} // projectionStep


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


/*!
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
  if (dfTilde != PETSC_NULL) {ierr = VecDestroy(&dfTilde); CHKERRQ(ierr);}
  if (tmp != PETSC_NULL) {ierr = VecDestroy(&tmp); CHKERRQ(ierr);}
  if (dlambda != PETSC_NULL) {ierr = VecDestroy(&dlambda); CHKERRQ(ierr);}
  if (rhs1_n != PETSC_NULL) {ierr = VecDestroy(&rhs1_n); CHKERRQ(ierr);}
  if (E != PETSC_NULL) {ierr = MatDestroy(&E); CHKERRQ(ierr);}
  if (ET != PETSC_NULL) {ierr = MatDestroy(&ET); CHKERRQ(ierr);}
  if (EBNET != PETSC_NULL) {ierr = MatDestroy(&EBNET); CHKERRQ(ierr);}
  if (G != PETSC_NULL) {ierr = MatDestroy(&G); CHKERRQ(ierr);}
  
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
#include "inline/calculateForces2.inl"
#include "inline/io.inl"


// dimensions specialization
template class LiEtAlSolver<2>;
template class LiEtAlSolver<3>;
