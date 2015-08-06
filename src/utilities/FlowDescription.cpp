/***************************************************************************//**
 * \file FlowDescription.cpp
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the methods of the class \c FlowDescription.
 */


#include "FlowDescription.h"

#include <iostream>
#include <fstream>

#include "yaml-cpp/yaml.h"


/**
 * \brief Converts \c std::string to \c Boundary.
 */
Boundary boundaryFromString(std::string s)
{
  if (s == "xMinus") return XMINUS;
  if (s == "xPlus")  return XPLUS;
  if (s == "yMinus") return YMINUS;
  if (s == "yPlus")  return YPLUS;
  if (s == "zMinus") return ZMINUS;
  if (s == "zPlus")  return ZPLUS;
  
  std::cout << "ERROR: Invalid boundary location!\n";
  exit(0);
} // boundaryFromString

/**
 * \brief Converts \c std::string to \c BCType.
 */
BCType bcTypeFromString(std::string s)
{
  if (s == "DIRICHLET") return DIRICHLET;
  if (s == "NEUMANN") return NEUMANN;
  if (s == "CONVECTIVE") return CONVECTIVE;
  if (s == "PERIODIC") return PERIODIC;
  
  std::cout << "ERROR: Invalid boundary condition type!\n";
  exit(0);
} // bcTypeFromString

FlowDescription::FlowDescription()
{
} // FlowDescription

/**
 * \brief Constructor -- Parses the input file.
 */
FlowDescription::FlowDescription(std::string fileName)
{
  initialize(fileName);
} // FlowDescription

/**
 * \brief Parses the input file and stores information about the flow.
 *
 * The file is parsed using YAML format.
 *
 * \param fileName path of the file to parse.
 */
void FlowDescription::initialize(std::string fileName)
{
  PetscInt    rank;
  
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank); // get rank of current process
  
  if (rank == 0) // read the input file only on process 0
  {
    YAML::Node nodes = YAML::LoadFile(fileName);
    const YAML::Node &node = nodes[0];

    dimensions = node["dimensions"].as<PetscInt>();
    nu = node["nu"].as<PetscReal>();

    initialVelocity[0] = node["initialVelocity"][0].as<PetscReal>();
    initialVelocity[1] = node["initialVelocity"][1].as<PetscReal>();
    initialVelocity[2] = node["initialVelocity"][2].as<PetscReal>(0.0);
    initialCustomField = (node["initialCustomField"].as<bool>(false)) ? PETSC_TRUE : PETSC_FALSE;

    perturbationAmplitude = node["initialPerturbation"][0].as<PetscReal>(0.0);
    perturbationFrequency = node["initialPerturbation"][1].as<PetscReal>(0);

    const YAML::Node &bcs = node["boundaryConditions"];
    Boundary location;
    // loop over the boundaries
    for (unsigned int i=0; i<bcs.size(); i++)
    {
      location = boundaryFromString(bcs[i]["location"].as<std::string>());
      bc[0][location].type = bcTypeFromString(bcs[i]["u"][0].as<std::string>());
      bc[0][location].value = bcs[i]["u"][1].as<PetscReal>();
      bc[1][location].type = bcTypeFromString(bcs[i]["v"][0].as<std::string>());
      bc[1][location].value = bcs[i]["v"][1].as<PetscReal>();
      // 2D cases: periodic boundary condition in the third direction
      bc[2][location].type = bcTypeFromString(bcs[i]["w"][0].as<std::string>("PERIODIC"));
      bc[2][location].value = bcs[i]["w"][1].as<PetscReal>(0.0);
    }
    
    // run some sanity checks on the input data
    PetscBool error = PETSC_FALSE;
    // if the boundary condition on one face is periodic,
    // then it should be periodic on the opposite face too
    // check u on the X and Y faces
    if (bc[0][XMINUS].type == PERIODIC && bc[0][XPLUS].type != PERIODIC) error = PETSC_TRUE;
    if (bc[0][XMINUS].type != PERIODIC && bc[0][XPLUS].type == PERIODIC) error = PETSC_TRUE;
    if (bc[0][YMINUS].type == PERIODIC && bc[0][YPLUS].type != PERIODIC) error = PETSC_TRUE;
    if (bc[0][YMINUS].type != PERIODIC && bc[0][YPLUS].type == PERIODIC) error = PETSC_TRUE;
    if (dimensions == 3)
    {
      if (bc[0][ZMINUS].type == PERIODIC && bc[0][ZPLUS].type != PERIODIC) error = PETSC_TRUE;
      if (bc[0][ZMINUS].type != PERIODIC && bc[0][ZPLUS].type == PERIODIC) error = PETSC_TRUE;
    }
    // if the boundary condition for one component of velocity is periodic,
    // then it should be periodic for the other components too
    // compare u and v on the X and Y faces
    if (bc[0][XMINUS].type == PERIODIC && bc[1][XMINUS].type != PERIODIC) error = PETSC_TRUE;
    if (bc[0][XPLUS].type == PERIODIC && bc[1][XPLUS].type != PERIODIC) error = PETSC_TRUE;
    if (bc[0][YMINUS].type == PERIODIC && bc[1][YMINUS].type != PERIODIC) error = PETSC_TRUE;
    if (bc[0][YPLUS].type == PERIODIC && bc[1][YPLUS].type != PERIODIC) error = PETSC_TRUE;
    if (dimensions == 3)
    {
      // compare u and w on the X and Y faces
      if (bc[0][XMINUS].type == PERIODIC && bc[2][XMINUS].type != PERIODIC) error = PETSC_TRUE;
      if (bc[0][XPLUS].type == PERIODIC && bc[2][XPLUS].type != PERIODIC) error = PETSC_TRUE;
      if (bc[0][YMINUS].type == PERIODIC && bc[2][YMINUS].type != PERIODIC) error = PETSC_TRUE;
      if (bc[0][YPLUS].type == PERIODIC && bc[2][YPLUS].type != PERIODIC) error = PETSC_TRUE;
      // compare u and v, u and w the on Z faces
      if (bc[0][ZMINUS].type == PERIODIC && bc[1][ZMINUS].type != PERIODIC) error = PETSC_TRUE;
      if (bc[0][ZMINUS].type == PERIODIC && bc[2][ZMINUS].type != PERIODIC) error = PETSC_TRUE;
      if (bc[0][ZPLUS].type == PERIODIC && bc[1][ZPLUS].type != PERIODIC) error = PETSC_TRUE;
      if (bc[0][ZPLUS].type == PERIODIC && bc[2][ZPLUS].type != PERIODIC) error = PETSC_TRUE;
    }
    if (error)
    {
      std::cout << "\nERROR: Check consistency of boundary conditions" << std::endl;
      exit(0);
    }
  }
  MPI_Barrier(PETSC_COMM_WORLD);
  
  // broadcast flow description to all processes
  MPI_Bcast(&dimensions, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
  MPI_Bcast(&nu, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
  MPI_Bcast(initialVelocity, 3, MPIU_REAL, 0, PETSC_COMM_WORLD);
  MPI_Bcast(&initialCustomField, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
  MPI_Bcast(&perturbationAmplitude, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
  MPI_Bcast(&perturbationFrequency, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
  
  // create custom MPI type to broadcast BC information
  MPI_Datatype bcInfoType, types[2];
  PetscInt blockcounts[2];
  MPI_Aint offsets[2];
  offsets[0] = offsetof(BoundaryCondition, type);
  types[0] = MPIU_INT;
  blockcounts[0] = 1;
  offsets[1] = offsetof(BoundaryCondition, value);
  types[1] = MPIU_REAL;
  blockcounts[1] = 1;
  MPI_Type_create_struct(2, blockcounts, offsets, types, &bcInfoType);
  MPI_Type_commit(&bcInfoType);
  MPI_Bcast(bc[0], 6, bcInfoType, 0, PETSC_COMM_WORLD);
  MPI_Bcast(bc[1], 6, bcInfoType, 0, PETSC_COMM_WORLD);
  MPI_Bcast(bc[2], 6, bcInfoType, 0, PETSC_COMM_WORLD);
  MPI_Type_free(&bcInfoType);
} // initialize