/***************************************************************************//**
 * \file CartesianMesh.h
 * \author Anush Krishnan (anus@bu.edu)
 * \brief Definition of the class \c CartesianMesh.
 */


#if !defined(CARTESIAN_MESH_H)
#define CARTESIAN_MESH_H

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

  // parse input file and initialize cartesian mesh
  CartesianMesh(std::string fileName);
  CartesianMesh();
  void initialize(std::string fileName);
};

#endif