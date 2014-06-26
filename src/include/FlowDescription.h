/***************************************************************************//**
* \file
* \brief Header to define class FlowDescription
*/

#if !defined(FLOW_DESCRIPTION_H)
#define FLOW_DESCRIPTION_H

#include "BoundaryCondition.h"
#include <petscsys.h>
#include <string>

/***************************************************************************//**
* \brief Stores information that describes the flow
*/
class FlowDescription
{
public:
	PetscInt          dimensions;             ///< number of dimensions 
	PetscReal         nu;                     ///< kinematic viscosity of the fluid
	PetscReal         initialVelocity[3];     ///< The initial velocity of the flow field
	PetscReal         initialPerturbation[3]; ///< The initial pertubation in the flow field
	BoundaryCondition bc[3][6];               ///< boundary conditions of the flow
	
	/***********************************************************************//**
	* \brief Reads the flow description from an input file
	*/
	FlowDescription(std::string fileName);
};

#endif
