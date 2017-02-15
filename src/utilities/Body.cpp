/*! Implementation of the methods of the class `Body`.
 * \file Body.cpp
 */


#include "Body.h"

#include <fstream>

#include "yaml-cpp/yaml.h"


/*!
 * \brief Constructor -- Reads the boundary coordinates from a given file.
 *
 * \param filePath Path of the file containing the boundary coordinates.
 */
template <PetscInt dim>
Body<dim>::Body(std::string filePath)
{
  readFromFile(filePath);
} // Body


/*!
 * \brief Reads the boundary coordinates from a given file.
 *
 * \param filePath Path of the file containing the boundary coordinates.
 */
template <PetscInt dim>
PetscErrorCode Body<dim>::readFromFile(std::string filePath)
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;

  ierr = PetscPrintf(PETSC_COMM_WORLD, "\nReading file %s... ", filePath.c_str()); CHKERRQ(ierr);

  std::ifstream infile(filePath.c_str());
  if (!infile.good())
  {
    ierr = PetscPrintf(PETSC_COMM_WORLD, "\nERROR: File '%s' does not exist\n", filePath.c_str()); CHKERRQ(ierr);
    exit(1);
  }
  infile >> numPoints;
  X.reserve(numPoints);
  Y.reserve(numPoints);
  PetscReal x, y;
  if (dim == 2)
  {
    for (PetscInt i=0; i<numPoints; i++)
    {
      infile >> x >> y;
      X.push_back(x);
      Y.push_back(y);
    }
  }
  else if (dim == 3)
  {
    PetscReal z;
    Z.reserve(numPoints);
    for (PetscInt i=0; i<numPoints; i++)
    {
      infile >> x >> y >> z;
      X.push_back(x);
      Y.push_back(y);
      Z.push_back(z);
    }
  }
  infile.close();

  ierr = PetscPrintf(PETSC_COMM_WORLD, "done.\n"); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // readFromFile


/*!
 * \brief Stores the indices of Eulerian cells owning a Lagrangian body points.
 * 
 * \param mesh Contains the information about the Cartesian grid
 */
template <PetscInt dim>
PetscErrorCode Body<dim>::registerCellOwners(CartesianMesh *mesh)
{
  I.reserve(numPoints);
  J.reserve(numPoints);
  if (dim == 3)
    K.reserve(numPoints);

  PetscInt i=0, j=0, k=0;

  // find cell owning first Lagrangian body point
  while (mesh->x[i+1] < X[0])
    i++;
  I.push_back(i);
  while (mesh->y[j+1] < Y[0])
    j++;
  J.push_back(j);
  if (dim == 3)
  {
    while (mesh->z[k+1] < Z[0])
      k++;
    K.push_back(k);
  }

  for (size_t l=1; l<X.size(); l++)
  {
    // x-component
    if (X[l] < X[l-1]) // body point located on left-side of previous one
    {
      while (mesh->x[i] > X[l])
        i--;
    }
    else // body point located on right-side of previous one
    {
      while (mesh->x[i+1] < X[l])
        i++;
    }
    I.push_back(i);
    // y-component
    if (Y[l] < Y[l-1]) // body point located on upper-side of previous one
    {
      while (mesh->y[j] > Y[l])
        j--;
    }
    else // body point located on lower-side of previous one
    {
      while (mesh->y[j+1] < Y[l])
        j++;
    }
    J.push_back(j);
    if (dim == 3)
    {
      // z-component
      if (Z[l] < Z[l-1]) // body point located on back-side of previous one
      {
        while (mesh->z[k] > Z[l])
          k--;
      }
      else // body point located on front-side of previous one
      {
        while (mesh->z[k+1] < Z[l])
          k++;
      }
      K.push_back(k);
    }
  }

  return 0;
} // registerCellOwners


/*!
 * \brief Registers the number of points on process.
 *
 * A process is represented by a box (xmin, xmax, ymin, ymax, zmin, zmax).
 *
 * \param box The box defining the process.
 */
template <PetscInt dim>
PetscErrorCode Body<dim>::registerNumPointsOnProcess(PetscReal (&box)[2*dim])
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  PetscMPIInt rank, size;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);

  numPointsOnProcess.resize(size);

  numPointsOnProcess[rank] = 0;
  for (PetscInt i=0; i<numPoints; i++)
  {
    if (box[0] <= X[i] && X[i] < box[1] &&
        box[2] <= Y[i] && Y[i] < box[3])
    {
      if (dim == 2 || dim == 3 && box[4] <= Z[i] && Z[i] < box[5])
      {
        numPointsOnProcess[rank]++;
      }
    }
  }

  ierr = MPI_Allgather(MPI_IN_PLACE, 1, MPIU_INT,
                       &numPointsOnProcess[0], 1, MPIU_INT,
                       PETSC_COMM_WORLD); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // registerNumPointsOnProcess


/*!
 * \brief Registers the number of points on the local process as well as the indices.
 *
 * A process is represented by a box (xmin, xmax, ymin, ymax, zmin, zmax).
 *
 * \param box The box defining the process.
 */
template <PetscInt dim>
PetscErrorCode Body<dim>::registerPointsOnProcess(PetscReal (&box)[2*dim])
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  idxPointsOnProcess.resize(numPoints);

  ierr = registerNumPointsOnProcess(box); CHKERRQ(ierr);

  PetscMPIInt rank, size;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);
  
  std::vector<PetscInt> offsets(size);
  offsets[rank] = 0;
  for (PetscMPIInt r=0; r<rank; r++)
  {
    offsets[rank] += numPointsOnProcess[r];
  }
  ierr = MPI_Allgather(MPI_IN_PLACE, 1, MPIU_INT,
                       &offsets[0], 1, MPIU_INT,
                       PETSC_COMM_WORLD); CHKERRQ(ierr);

  PetscInt index = offsets[rank];
  for (PetscInt i=0; i<numPoints; i++)
  {
    if (box[0] <= X[i] && X[i] < box[1] &&
        box[2] <= Y[i] && Y[i] < box[3])
    {
      if (dim == 2 || (dim == 3 && box[4] <= Z[i] && Z[i] < box[5]))
      {
        idxPointsOnProcess[index] = i;
        index++;
      }
    }
  }
  ierr = MPI_Allgatherv(MPI_IN_PLACE, numPointsOnProcess[rank], MPIU_INT,
                        &idxPointsOnProcess[0], &numPointsOnProcess[0], &offsets[0], MPIU_INT,
                        PETSC_COMM_WORLD); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // registerPointsOnProcess


/*!
 * \brief Registers the global index on points on local process.
 *
 * The global index represents the position of a Lagrangian force in the global
 * vector lambda.
 * \param globalIdx Offset value to increment from.
 */
template <PetscInt dim>
PetscErrorCode Body<dim>::registerGlobalIdxPoints(PetscInt &globalIdx)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  globalIdxPoints.resize(numPoints);

  PetscMPIInt rank, size;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);
  
  std::vector<PetscInt> offsets(size);
  offsets[rank] = 0;
  for (PetscMPIInt r=0; r<rank; r++)
  {
    offsets[rank] += numPointsOnProcess[r];
  }
  ierr = MPI_Allgather(MPI_IN_PLACE, 1, MPIU_INT,
                       &offsets[0], 1, MPIU_INT,
                       PETSC_COMM_WORLD); CHKERRQ(ierr);

  for (PetscInt i=offsets[rank]; i<offsets[rank]+numPointsOnProcess[rank]; i++)
  {
    globalIdxPoints[i] = globalIdx;
    globalIdx += dim;
  }

  ierr = MPI_Allgatherv(MPI_IN_PLACE, numPointsOnProcess[rank], MPIU_INT,
                        &globalIdxPoints[0], &numPointsOnProcess[0], &offsets[0], MPIU_INT,
                        PETSC_COMM_WORLD); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // registerGlobalIdxPoints


// dimensions specialization
template class Body<2>;
template class Body<3>;
