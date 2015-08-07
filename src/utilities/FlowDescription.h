/***************************************************************************//**
 * \file FlowDescription.h
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Definition of the class \c FlowDescription.
 */


#if !defined(FLOW_DESCRIPTION_H)
#define FLOW_DESCRIPTION_H

// #include "BoundaryCondition.h"
#include "types.h"

#include <string>

#include <petscsys.h>


/**
 * \class FlowDescription
 * \brief Stores information that describes the flow
 */
class FlowDescription
{
public:
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

  PetscInt dimensions;              ///< number of dimensions 
  PetscReal nu;                     ///< kinematic viscosity of the fluid
  PetscReal initialVelocity[3];     ///< initial velocity of the flow field
  PetscBool initialCustomField;     ///< flag to start from an initial custom field
  PetscReal perturbationAmplitude;  ///< amplitude of the Taylor-Green vortex perturbation
  PetscReal perturbationFrequency;  ///< frequency of the Taylor-Green vortex perturbation
  BoundaryCondition bc[3][6];       ///< boundary conditions of the flow
  
  // contructors
  FlowDescription();
  FlowDescription(std::string directory);
  // destructor
  ~FlowDescription();
  // parse input file and store description of the flow
  void initialize(std::string filePath);

}; // FlowDescription

#endif