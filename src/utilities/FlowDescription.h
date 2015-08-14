/***************************************************************************//**
 * \file FlowDescription.h
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Definition of the class \c FlowDescription.
 */


#if !defined(FLOW_DESCRIPTION_H)
#define FLOW_DESCRIPTION_H

#include "types.h"

#include <string>

#include <petscsys.h>


/**
 * \class FlowDescription
 * \brief Stores information that describes the flow
 */
template <PetscInt dim>
class FlowDescription
{
public:
  /**
   * \class Boundary
   * \brief Stores info about a boundary (type of boundary condition type and value).
   */
  class Boundary
  {
  public:
    BoundaryType type; ///< type of boundary condition
    PetscReal value;   ///< value at the boundary
  }; // Boundary

  PetscReal nu;                    ///< kinematic viscosity of the fluid
  PetscReal initialVelocity[dim];  ///< initial velocity of the flow field
  PetscBool initialCustomField;    ///< flag to start from an initial custom field
  PetscReal perturbationAmplitude; ///< amplitude of the Taylor-Green vortex perturbation
  PetscReal perturbationFrequency; ///< frequency of the Taylor-Green vortex perturbation
  Boundary boundaries[2*dim][dim]; ///< boundary conditions of the flow
  
  // contructors
  FlowDescription();
  FlowDescription(std::string directory);
  // destructor
  ~FlowDescription();
  // parse input file and store description of the flow
  void initialize(std::string filePath);
  // run sanity checks about periodic boundary conditions
  void checkPeriodicity();
  // print description of the flow
  PetscErrorCode printInfo();

}; // FlowDescription

#endif