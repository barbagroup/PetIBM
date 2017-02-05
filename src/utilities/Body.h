/*! Definition of the class `Body`.
 * \file Body.h
 */


#if !defined(BODY_H)
#define BODY_H

#include "CartesianMesh.h"

#include <vector>
#include <string>

#include <petscsys.h>


/**
 * \class Body.h
 * \brief Contains information about an immersed boundary.
 */
template <PetscInt dim>
class Body
{
public:
  PetscInt numPoints; ///< number of points constituting the body
  std::vector<PetscReal> X, ///< x-coordinate of all body points
                         Y, ///< y-coordinate of all body points
                         Z; ///< z-coordinate of all body points
  std::vector<PetscInt> I, ///< x-index of Eulerian cells owning the body points
                        J, ///< x-index of Eulerian cells owning the body points
                        K; ///< x-index of Eulerian cells owning the body points
  PetscReal forces[dim]; ///< Force vector acting on the body
  std::vector<PetscInt> numPointsOnProcess; ///< number of body points on each process
  std::vector<PetscInt> idxPointsOnProcess; ///< index of body points per process
  std::vector<PetscInt> globalIdxPoints; ///< local-to-global mapping

  // constructors
  Body(){ };
  Body(std::string filePath);
  // destructor
  ~Body(){ };

  // read the body coordinates from file
  PetscErrorCode readFromFile(std::string filePath);
  // register the indices of cells owning a Lagrangian body point
  PetscErrorCode registerCellOwners(CartesianMesh *mesh);
  // register the number of body points in a given box
  PetscErrorCode registerNumPointsOnProcess(PetscReal (&box)[2*dim]);
  // register the number and index of body points in a given box
  PetscErrorCode registerPointsOnProcess(PetscReal (&box)[2*dim]);
  // register the local-to-global mapping
  PetscErrorCode registerGlobalIdxPoints(PetscInt &offset);

}; // Body

#endif
