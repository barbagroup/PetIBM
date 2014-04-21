#if !defined(BOUNDARY_CONDITION_H)
#define BOUNDARY_CONDITION_H

#include <petscsys.h>
#include "types.h"

class BoundaryCondition
{
public:
	BCType     type; ///< Type of boundary condition
	PetscReal  value; ///< Numerical value associated with the boundary condition
};

#endif