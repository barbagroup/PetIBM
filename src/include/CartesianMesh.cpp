/***************************************************************************//**
 * \file CartesianMesh.cpp
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the class \c CartesianMesh.
 */


#include "CartesianMesh.h"

#include <fstream>

#include "yaml-cpp/yaml.h"


CartesianMesh::CartesianMesh()
{
}

/**
 * \brief Constructor -- Parses input file.
 */
CartesianMesh::CartesianMesh(std::string fileName)
{
  initialize(fileName);
}

/**
 * \brief Parses the input file using YAML format and discretizes the domain.
 *
 * \param fileName path of the file containing the mesh parameters
 */
void CartesianMesh::initialize(std::string fileName)
{
  PetscInt rank;
  
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank); // get rank of the current process
  
  nx = 0;
  ny = 0;
  nz = 0;
  
  // first pass of the input file to get number of cells in each direction
  if (rank == 0)
  {
    YAML::Node nodes = YAML::LoadFile(fileName);
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
        nx = numCells;
      else if (directio == "y")
        ny = numCells;
      else if (direction == "z")
        nx += numCells;
    }

  MPI_Barrier(PETSC_COMM_WORLD);
  
  // broadcast number of cells to all processes
  MPI_Bcast(&nx, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
  MPI_Bcast(&ny, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
  MPI_Bcast(&nz, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
  
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
  if (rank == 0)
  {
    PetscInt numCells,
             first;
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
          for (unsigned int j=first; j<first+numCells; j++)
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
          for (unsigned int j=first; j<first+numCells; j++)
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
  }
  MPI_Barrier(PETSC_COMM_WORLD);
  
  // broadcast vectors to all processes
  MPI_Bcast(&x.front(), nx+1, MPIU_REAL, 0, PETSC_COMM_WORLD);
  MPI_Bcast(&dx.front(), nx, MPIU_REAL, 0, PETSC_COMM_WORLD);
  MPI_Bcast(&y.front(), ny+1, MPIU_REAL, 0, PETSC_COMM_WORLD);
  MPI_Bcast(&dy.front(), ny, MPIU_REAL, 0, PETSC_COMM_WORLD);
  if (nz > 0)
  {
    MPI_Bcast(&z.front(), nz+1, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(&dz.front(), nz, MPIU_REAL, 0, PETSC_COMM_WORLD);
  }
}