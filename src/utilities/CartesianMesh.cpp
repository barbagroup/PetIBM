/***************************************************************************//**
 * \file CartesianMesh.cpp
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the class \c CartesianMesh.
 */


#include "CartesianMesh.h"

#include <fstream>

#include "yaml-cpp/yaml.h"


/**
 * \brief Constructor.
 */
CartesianMesh::CartesianMesh()
{
} // CartesianMesh


/**
 * \brief Constructor -- Parses input file.
 */
CartesianMesh::CartesianMesh(std::string directory)
{
  initialize(directory + "/cartesianMesh.yaml");
} // CartesianMesh


/**
 * \brief Destructor
 */
CartesianMesh::~CartesianMesh()
{
} // ~CartesianMesh


/**
 * \brief Parses the input file using YAML format and discretizes the domain.
 *
 * \param filePath path of the file containing the mesh parameters
 */
void CartesianMesh::initialize(std::string filePath)
{
  PetscPrintf(PETSC_COMM_WORLD, "\nParsing file %s... ", filePath.c_str());
  
  nx = 0;
  ny = 0;
  nz = 0;
  
  // first pass of the input file to get number of cells in each direction
  YAML::Node nodes = YAML::LoadFile(filePath);
  PetscInt numCells;
  std::string direction;
  for (unsigned int i=0; i<nodes.size(); i++)
  {
    numCells = 0;
    direction = nodes[i]["direction"].as<std::string>();
    const YAML::Node &subDomains = nodes[i]["subDomains"];
    for (unsigned int j=0; j<subDomains.size(); j++)
      numCells += subDomains[j]["cells"].as<PetscInt>();
    if (direction == "x")
      nx += numCells;
    else if (direction == "y")
      ny += numCells;
    else if (direction == "z")
      nz += numCells;
  }
  
  // allocate memory
  // number of nodes will be 1 greater than the number of cells
  x.resize(nx+1);
  dx.resize(nx);
  y.resize(ny+1);
  dy.resize(ny);
  if (nz > 0)
  {
    z.resize(nz+1);
    dz.resize(nz);
  }
  
  // second pass of the input file
  // to calculate the coordinates of the nodes, and the cell widths
  PetscInt first;
  PetscReal start,
            end,
            stretchRatio,
            h;

  // loop over each direction
  for (unsigned int k=0; k<nodes.size(); k++)
  {
    direction = nodes[k]["direction"].as<std::string>();
    start = nodes[k]["start"].as<PetscReal>();

    const YAML::Node &subDomains = nodes[k]["subDomains"];
    first = 0;
    if (direction == "x") x[first] = start;
    else if (direction == "y") y[first] = start;
    else if (direction == "z") z[first] = start;
    for (unsigned int i=0; i<subDomains.size(); i++)
    {
      end = subDomains[i]["end"].as<PetscReal>();
      numCells = subDomains[i]["cells"].as<PetscInt>();
      stretchRatio = subDomains[i]["stretchRatio"].as<PetscReal>();

      if (fabs(stretchRatio-1.0) < 1.0E-06) // uniform discretization
      {
        h = (end - start)/numCells;
        for (int j=first; j<first+numCells; j++)
        {
          if (direction == "x")
          {
            dx[j] = h;
            x[j+1] = x[j] + dx[j];
          }
          else if (direction == "y")
          {
            dy[j] = h;
            y[j+1] = y[j] + dy[j];
          }
          else if (direction == "z")
          {
            dz[j] = h;
            z[j+1] = z[j] + dz[j];
          }
        }
      }
      else // stretched discretization
      {
        // width of the first cell in the subdomain
        h = (end - start)*(stretchRatio-1)/(pow(stretchRatio, numCells)-1);
        for (int j=first; j<first+numCells; j++)
        {
          // obtain the widths of subsequent cells by
          // multiplying by the stretching ratio
          if (direction == "x")
          {
            dx[j]  = h*pow(stretchRatio, j-first);
            x[j+1] = x[j] + dx[j];
          }
          else if (direction == "y")
          {
            dy[j]  = h*pow(stretchRatio, j-first);
            y[j+1] = y[j] + dy[j];
          }
          else if (direction == "z")
          {
            dz[j]  = h*pow(stretchRatio, j-first);
            z[j+1] = z[j] + dz[j];
          }
        }
      }
      // the first node in the next subdomain is the last node in the current subdomain
      first += numCells;
      start = end;
    }
  }

  PetscPrintf(PETSC_COMM_WORLD, "done.\n");

} // initialize


/**
 * \brief Prints information about the Cartesian mesh.
 */
PetscErrorCode CartesianMesh::printInfo()
{
  PetscErrorCode ierr;
  ierr = PetscPrintf(PETSC_COMM_WORLD, "\n---------------------------------------\n"); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Cartesian grid\n"); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "---------------------------------------\n"); CHKERRQ(ierr);
  if (nz > 0)
  {
    ierr = PetscPrintf(PETSC_COMM_WORLD, "number of cells: %d x %d x %d\n", nx, ny, nz); CHKERRQ(ierr);
  }
  else
  {
    ierr = PetscPrintf(PETSC_COMM_WORLD, "number of cells: %d x %d\n", nx, ny); CHKERRQ(ierr);
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD, "---------------------------------------\n"); CHKERRQ(ierr);

  return 0;
} // printInfo