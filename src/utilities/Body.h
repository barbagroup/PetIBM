/***************************************************************************//**
 * \file Body.h
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \brief Definition of the class `Body`.
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

  // parse input file and store coordinates of immersed body
  PetscErrorCode initialize(std::string filePath);
  // get indices of cells owning a Lagrangian body point
  PetscErrorCode getCellOwners(CartesianMesh *mesh);
  // print info related to body
  PetscErrorCode printInfo();

}; // Body

#endif