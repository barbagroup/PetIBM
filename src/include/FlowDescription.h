#if !defined(FLOW_DESCRIPTION_H)
#define FLOW_DESCRIPTION_H

#include "BoundaryCondition.h"
#include <petscsys.h>
#include <string>

class FlowDescription
{
public:
	PetscInt  dimensions;
	PetscReal nu;
	PetscReal initialVelocity[3];
	BoundaryCondition bc[3][6];
	
	FlowDescription(std::string fileName);
};

#endif
