/***************************************************************************//**
 * \file Body.cpp
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \brief Implementation of the methods of the class `Body`.
 */


#include "Body.h"

#include <fstream>

#include "yaml-cpp/yaml.h"


/**
 * \brief Constructor.
 */
template <PetscInt dim>
Body<dim>::Body()
{
} // Body


/**
 * \brief Constructor -- Parses the input file `bodies.yaml`.
 *
 * \param directory Directory of the simulation
 */
template <PetscInt dim>
Body<dim>::Body(std::string directory)
{
  initialize(directory + "/bodies.yaml");
} // Body


/**
 * \brief Destructor.
 */
template <PetscInt dim>
Body<dim>::~Body()
{
} // ~Body


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
      exit(0);
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

  for (PetscInt l=1; l<X.size(); l++)
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


// dimensions specialization
template class Body<2>;
template class Body<3>;