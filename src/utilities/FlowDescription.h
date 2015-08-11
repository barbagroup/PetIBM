/***************************************************************************//**
 * \file FlowDescription.h
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Definition of the class \c FlowDescription.
 */


#if !defined(FLOW_DESCRIPTION_H)
#define FLOW_DESCRIPTION_H

#include "BoundaryCondition.h"

#include <string>

#include <petscsys.h>


/**
 * \class FlowDescription
 * \brief Stores information that describes the flow
 */
class FlowDescription
{
public:
  PetscInt dimensions;              ///< number of dimensions 
  PetscReal nu;                     ///< kinematic viscosity of the fluid
  PetscReal initialVelocity[3];     ///< initial velocity of the flow field
  PetscBool initialCustomField;     ///< flag to start from an initial custom field
  PetscReal perturbationAmplitude;  ///< amplitude of the Taylor-Green vortex perturbation
  PetscReal perturbationFrequency;  ///< frequency of the Taylor-Green vortex perturbation
  BoundaryCondition bc[3][6];       ///< boundary conditions of the flow
  
  // parse the input file and store information about the flow
  FlowDescription(std::string fileName);
  FlowDescription();
  void initialize(std::string fileName);

}; // FlowDescription

#endif