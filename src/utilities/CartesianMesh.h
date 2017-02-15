/***************************************************************************//**
 * \file CartesianMesh.h
 * \author Anush Krishnan (anus@bu.edu)
 * \brief Definition of the class `CartesianMesh`.
 */


#if !defined(CARTESIAN_MESH_H)
#define CARTESIAN_MESH_H

#include "types.h"

#include <string>
#include <vector>

#include <petscsys.h>


/**
 * \brief Stores information about the cartesian mesh.
 */
class CartesianMesh
{
public:
  PetscInt nx, ///< number of cells in the x-direction
           ny, ///< number of cells in the y-direction
           nz; ///< number of cells in the z-direction
  
  std::vector<PetscReal> x, ///< x-coordinates of the nodes
                         y, ///< y-coordinates of the nodes
                         z; ///< z-coordinates of the nodes
  
  std::vector<PetscReal> dx, ///< cell-widths along the x-direction
                         dy, ///< cell-widths along the y-direction 
                         dz; ///< cell-widths along the z-direction

  // constructors
  CartesianMesh();
  CartesianMesh(std::string filePath);
  // destructor
  ~CartesianMesh();
  // parse input file and create Cartesian mesh
  void initialize(std::string filePath);
  // write grid points into file
  PetscErrorCode write(std::string filePath);
#ifdef PETSC_HAVE_HDF5
  PetscErrorCode write(std::string filePath, StaggeredMode mode, BoundaryType type);
#endif
  // print information about Cartesian mesh
  PetscErrorCode printInfo();

}; // CartesianMesh

#endif