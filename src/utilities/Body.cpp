/***************************************************************************//**
 * \file Body.cpp
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \brief Implementation of the methods of the class `Body`.
 */


#include "Body.h"

#include <fstream>

#include "yaml-cpp/yaml.h"


/**
 * \brief Constructor -- Parses the YAMLinput file containing bodies information.
 *
 * \param filePath Path of the file to parse with YAML-CPP
 */
template <PetscInt dim>
Body<dim>::Body(std::string filePath)
{
  char path[PETSC_MAX_PATH_LEN];
  PetscBool found;
  PetscOptionsGetString(NULL, NULL, "-bodies", path, sizeof(path), &found);
  if (found)
    filePath = std::string(path);
  initialize(filePath);
} // Body

/*!
 * \brief Constructor -- Reads the boundary coordinates from a given file.
 *
 * \param filePath Path of the file containing the boundary coordinates.
 */
// template <PetscInt dim>
// Body<dim>::Body(std::string filePath)
// {
//   readFromFile(filePath);
// } // Body


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
  PetscInt x, y;
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
    PetscInt z;
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
 * \brief Parses the input file and stores information about the flow.
 *
 * The file is parsed using YAML format.
 *
 * \param filePath Path of the file to parse
 */
template <PetscInt dim>
PetscErrorCode Body<dim>::initialize(std::string filePath)
{
  PetscErrorCode ierr;

  ierr = PetscPrintf(PETSC_COMM_WORLD, "\nParsing file %s... ", filePath.c_str()); CHKERRQ(ierr);

  YAML::Node nodes = YAML::LoadFile(filePath);
  const YAML::Node &node = nodes[0];

  std::string type = node["type"].as<std::string>();

  if (type == "points")
  {
    PetscReal x, y, z;
    size_t last = filePath.find_last_of("/");
    std::string directory = filePath.substr(0, last+1);
    std::string pointsFilePath = directory + node["pointsFile"].as<std::string>();
    std::ifstream infile(pointsFilePath.c_str());
    if (!infile.good())
    {
      ierr = PetscPrintf(PETSC_COMM_WORLD, "\nERROR: File '%s' does not exist\n", pointsFilePath.c_str()); CHKERRQ(ierr);
      exit(1);
    }
    infile >> numPoints;
    X.reserve(numPoints);
    Y.reserve(numPoints);
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
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD, "done.\n"); CHKERRQ(ierr);

  return 0;
} // initialize


/**
 * \brief Gets the indices of cells owning a Lagrangian body points.
 * 
 * \param mesh Contains the information about the Cartesian grid
 */
template <PetscInt dim>
PetscErrorCode Body<dim>::getCellOwners(CartesianMesh *mesh)
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
} // getCellOwners


/*!
 * \brief Counts the number of points in a given box.
 *
 * \param box Box defined by (xmin, xmax, ymin, ymax, zmin, zmax).
 *
 * \returns The number of points in the box.
 */
template <PetscInt dim>
PetscErrorCode Body<dim>::countNumberPointsInBox(PetscReal (&box)[2*dim], PetscInt &n)
{
  return 0;
} // countNumberPointsInBox

// two-dimensional specialization
template <>
PetscErrorCode Body<2>::countNumberPointsInBox(PetscReal (&box)[4], PetscInt &n)
{
  PetscFunctionBeginUser;

  for (PetscInt i=0; i<numPoints; i++)
  {
    if (box[0] <= X[i] && X[i] < box[1] &&
        box[2] <= Y[i] && Y[i] < box[3])
      n++;
  }

  PetscFunctionReturn(0);
} // countNumberPointsInBox


// three-dimensional specialization
template <>
PetscErrorCode Body<3>::countNumberPointsInBox(PetscReal (&box)[6], PetscInt &n)
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
} // countNumberPointsInBox


/*!
 * \brief Gets the index of points inside a given box.
 *
 * \param box Box defined by (xmin, xmax, ymin, ymax, zmax, zmin).
 * \param indices The vector of indices to fill.
 */
template <PetscInt dim>
PetscErrorCode Body<dim>::getIndexPointsInBox(PetscReal (&box)[2*dim], std::vector<PetscInt> indices)
{
  return 0;
} // getIndexPointsInBox


// two-dimensional specialization
template <>
PetscErrorCode Body<2>::getIndexPointsInBox(PetscReal (&box)[4], std::vector<PetscInt> indices)
{
  PetscFunctionBeginUser;

  for (PetscInt i=0; i<numPoints; i++)
  {
    if (box[0] <= X[i] && X[i] < box[1] &&
        box[2] <= Y[i] && Y[i] < box[3])
      indices.push_back(i);
  }

  PetscFunctionReturn(0);
} // getIndexPointsInBox


// three-dimensional specialization
template <>
PetscErrorCode Body<3>::getIndexPointsInBox(PetscReal (&box)[6], std::vector<PetscInt> indices)
{
  PetscFunctionBeginUser;

  for (PetscInt i=0; i<numPoints; i++)
  {
    if (box[0] <= X[i] && X[i] < box[1] &&
        box[2] <= Y[i] && Y[i] < box[3] &&
        box[4] <= Z[i] && Z[i] < box[5])
      indices.push_back(i);
  }

  PetscFunctionReturn(0);
} // getIndexPointsInBox


// dimensions specialization
template class Body<2>;
template class Body<3>;
