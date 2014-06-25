/***************************************************************************//**
* \file
* \brief Header to define class BoundaryCondition
*/

#if !defined(BOUNDARY_CONDITION_H)
#define BOUNDARY_CONDITION_H

#include <petscsys.h>
#include "types.h"

/***************************************************************************//**
* \class BoundaryCondition
* \brief Store the type of boundary condition and the associated value
*/
class BoundaryCondition
{
public:
	BCType     type;  ///< Type of boundary condition
	PetscReal  value; ///< Numerical value associated with the boundary condition
};

#endif