/***************************************************************************//**
 * \file Body.cpp
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \brief Implementation of the methods of the class `Body`.
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
 * \brief Resizes the vector with the number of processors.
 */
template <PetscInt dim>
PetscErrorCode Body<dim>::initialize()
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  PetscMPIInt numProcs;
  ierr = MPI_Comm_size(PETSC_COMM_WORLD, &numProcs); CHKERRQ(ierr);

  localNumPoints.resize(numProcs);
  localIndexPoints.resize(numProcs);

  PetscFunctionReturn(0);
} // initialize


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


/**
 * \brief Stores the indices of Eulerian cells owning a Lagrangian body points.
 * 
 * \param mesh Contains the information about the Cartesian grid
 */
template <PetscInt dim>
PetscErrorCode Body<dim>::setCellOwners(CartesianMesh *mesh)
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
} // setCellOwners


/*!
 * \brief Gets the number of points in a given box.
 *
 * \param box Box defined by (xmin, xmax, ymin, ymax, zmin, zmax).
 * \param n The number of points in the box.
 */
template <PetscInt dim>
PetscErrorCode Body<dim>::getNumPointsInBox(PetscReal (&box)[2*dim], PetscInt &n)
{
  return 0;
} // getNumPointsInBox

// two-dimensional specialization
template <>
PetscErrorCode Body<2>::getNumPointsInBox(PetscReal (&box)[4], PetscInt &n)
{
  PetscFunctionBeginUser;

  for (PetscInt i=0; i<numPoints; i++)
  {
    if (box[0] <= X[i] && X[i] < box[1] &&
        box[2] <= Y[i] && Y[i] < box[3])
      n++;
  }

  PetscFunctionReturn(0);
} // getNumPointsInBox


// three-dimensional specialization
template <>
PetscErrorCode Body<3>::getNumPointsInBox(PetscReal (&box)[6], PetscInt &n)
{
  PetscFunctionBeginUser;

  for (PetscInt i=0; i<numPoints; i++)
  {
    if (box[0] <= X[i] && X[i] < box[1] &&
        box[2] <= Y[i] && Y[i] < box[3] &&
        box[4] <= Z[i] && Z[i] < box[5])
      n++;
  }

  PetscFunctionReturn(0);
} // getNumPointsInBox


/*!
 * \brief Stores the index of boundary points that on a given process
          (defined by a box).
 *
 * \param procIdx The process id.
 * \param box The box defining the process.
 */
template <PetscInt dim>
PetscErrorCode Body<dim>::setLocalIndexPoints(PetscInt procIdx, PetscReal (&box)[2*dim])
{
  return 0;
} // setLocalIndexPoints


// two-dimensional specialization
template <>
PetscErrorCode Body<2>::setLocalIndexPoints(PetscInt procIdx, PetscReal (&box)[4])
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  ierr = getNumPointsInBox(box, localNumPoints[procIdx]); CHKERRQ(ierr);
  localIndexPoints[procIdx].reserve(localNumPoints[procIdx]);
  for (PetscInt i=0; i<numPoints; i++)
  {
    if (box[0] <= X[i] && X[i] < box[1] &&
        box[2] <= Y[i] && Y[i] < box[3])
      localIndexPoints[procIdx].push_back(i);
  }

  PetscFunctionReturn(0);
} // setLocalIndexPoints


// three-dimensional specialization
template <>
PetscErrorCode Body<3>::setLocalIndexPoints(PetscInt procIdx, PetscReal (&box)[6])
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  ierr = getNumPointsInBox(box, localNumPoints[procIdx]); CHKERRQ(ierr);
  localIndexPoints[procIdx].reserve(localNumPoints[procIdx]);
  for (PetscInt i=0; i<numPoints; i++)
  {
    if (box[0] <= X[i] && X[i] < box[1] &&
        box[2] <= Y[i] && Y[i] < box[3] &&
        box[4] <= Z[i] && Z[i] < box[5])
      localIndexPoints[procIdx].push_back(i);
  }
  
  PetscFunctionReturn(0);
} // setLocalIndexPoints


// dimensions specialization
template class Body<2>;
template class Body<3>;
