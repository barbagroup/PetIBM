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
  std::vector<PetscInt> localNumPoints;
  std::vector<std::vector<PetscInt> > localIndexPoints;
  std::vector<PetscReal> X, Y, Z;
  std::vector<PetscInt> I, J, K;
  PetscReal forces[dim];

  // constructors
  Body(){ };
  Body(std::string filePath);
  // destructor
  ~Body(){ };

  // initialization
  PetscErrorCode initialize();
  // read boundary coordinates from give file
  PetscErrorCode readFromFile(std::string filePath);
  // store indices of cells owning a Lagrangian body point
  PetscErrorCode setCellOwners(CartesianMesh *mesh);
  // get number of points in a given box
  PetscErrorCode getNumPointsInBox(PetscReal (&box)[2*dim], PetscInt &n);
  // store index of points for each process
  PetscErrorCode setLocalIndexPoints(PetscInt procIdx, PetscReal (&box)[2*dim]);
  // print info related to body
  PetscErrorCode printInfo();

}; // Body

#endif
