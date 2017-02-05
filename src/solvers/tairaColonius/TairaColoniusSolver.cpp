/*! Implementation of the methods of the class `TairaColoniusSolver`.
 * \file TairaColoniusSolver.cpp
 */


#include "TairaColoniusSolver.h"

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
TairaColoniusSolver<dim>::TairaColoniusSolver(CartesianMesh *cartesianMesh, 
                                              FlowDescription<dim> *flowDescription, 
                                              SimulationParameters *simulationParameters)
                        : NavierStokesSolver<dim>::NavierStokesSolver(cartesianMesh, 
                                                                      flowDescription, 
                                                                      simulationParameters)
{
  bda = PETSC_NULL;
  nullSpaceVec = PETSC_NULL;
} // TairaColoniusSolver


/*!
 * \brief Initializes the solver.
 */
template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::initialize()
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  ierr = PetscLogStagePush(NavierStokesSolver<dim>::stageInitialize); CHKERRQ(ierr);

  ierr = initializeBodies(); CHKERRQ(ierr);
  ierr = createDMs(); CHKERRQ(ierr);
  ierr = createVecs(); CHKERRQ(ierr);
  ierr = createGlobalMappingBodies(); CHKERRQ(ierr);

  ierr = NavierStokesSolver<dim>::initializeCommon(); CHKERRQ(ierr);

  ierr = PetscLogStagePop(); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // initialize


/*!
 * \brief Gets the total number of Lagrangian points.
 *
 * \param n The integer to return.
 */
template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::getNumLagPoints(PetscInt &n)
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
PetscErrorCode TairaColoniusSolver<dim>::getNumLagPointsOnProcess(std::vector<PetscInt> &numOnProcess)
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

  ierr = MPI_Allgather(&numOnProcess[rank], 1, MPIU_INT,
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
PetscErrorCode TairaColoniusSolver<dim>::registerLagPointsOnProcess()
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
PetscErrorCode TairaColoniusSolver<dim>::finalize()
{
  PetscErrorCode ierr;

  ierr = NavierStokesSolver<dim>::finalize();
  // DMs
  if (bda != PETSC_NULL) {ierr = DMDestroy(&bda); CHKERRQ(ierr);}
  // Vecs
  if (nullSpaceVec != PETSC_NULL)    {ierr = VecDestroy(&nullSpaceVec); CHKERRQ(ierr);}

  return 0;
}  // finalize


#include "inline/initializeBodies.inl"
#include "inline/createDMs.inl"
#include "inline/createVecs.inl"
#include "inline/generateBNQ.inl"
#include "inline/generateR2.inl"
#include "inline/createGlobalMappingBodies.inl"
#include "inline/setNullSpace.inl"
#include "inline/calculateForces.inl"
#include "inline/io.inl"


// dimensions specialization
template class TairaColoniusSolver<2>;
template class TairaColoniusSolver<3>;
