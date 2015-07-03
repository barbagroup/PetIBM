/***************************************************************************//**
 * \file BoundaryCondition.h
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Definition of the class \c BoundaryCondition.
 */


#if !defined(BOUNDARY_CONDITION_H)
#define BOUNDARY_CONDITION_H

#include "types.h"

#include <petscsys.h>


/**
 * \class BoundaryCondition
 * \brief Stores the type of boundary condition and its associated value.
 */
class BoundaryCondition
{
public:
  BCType type;     ///< Type of boundary condition
  PetscReal value; ///< Numerical value associated with the boundary condition

}; // BoundaryCondition

#endif