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
  PetscInt numPoints;
  std::vector<PetscReal> X, Y, Z;
  std::vector<PetscInt> I, J, K;

  // constructors
  Body(){ };
  Body(std::string filePath);
  // destructor
  ~Body(){ };

  // read boundary coordinates from give file
  PetscErrorCode readFromFile(std::string filePath);
  // get indices of cells owning a Lagrangian body point
  PetscErrorCode getCellOwners(CartesianMesh *mesh);
  // count number of points in a given box
  PetscErrorCode countNumberPointsInBox(PetscReal (&box)[2*dim], PetscInt &n);
  // get index of points inside a given box
  PetscErrorCode getIndexPointsInBox(PetscReal (&box)[2*dim], std::vector<PetscInt> &indices);
  // print info related to body
  PetscErrorCode printInfo();

}; // Body

#endif