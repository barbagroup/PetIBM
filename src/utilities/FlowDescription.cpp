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
 * \brief Constructor.
 */
template <PetscInt dim>
FlowDescription<dim>::FlowDescription()
{
} // FlowDescription


/**
 * \brief Constructor -- Parses the input file.
 */
template <PetscInt dim>
FlowDescription<dim>::FlowDescription(std::string directory)
{
  initialize(directory + "/flowDescription.yaml");
} // FlowDescription


/**
 * \brief Destructor
 */
template <PetscInt dim>
FlowDescription<dim>::~FlowDescription()
{
} // ~FlowDescription


/**
 * \brief Parses the input file and stores information about the flow.
 *
 * The file is parsed using YAML format.
 *
 * \param filePath path of the file to parse.
 */
template <PetscInt dim>
void FlowDescription<dim>::initialize(std::string filePath)
{
  PetscPrintf(PETSC_COMM_WORLD, "\nParsing file %s... ", filePath.c_str());
  
  // get rank of current process
  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  
  if (rank == 0) // parse the input file only on process 0
  {
    YAML::Node nodes = YAML::LoadFile(filePath);
    const YAML::Node &node = nodes[0];

    nu = node["nu"].as<PetscReal>();

    for (PetscInt i=0; i<dim; i++)
    {
      initialVelocity[i] = node["initialVelocity"][i].as<PetscReal>();
    }

    initialCustomField = (node["initialCustomField"].as<bool>(false)) ? PETSC_TRUE : PETSC_FALSE;

    perturbationAmplitude = node["initialPerturbation"][0].as<PetscReal>(0.0);
    perturbationFrequency = node["initialPerturbation"][1].as<PetscReal>(0);

    const YAML::Node &bcs = node["boundaryConditions"];
    PetscBool checkDimensions = (bcs.size() == 2*dim) ? PETSC_TRUE : PETSC_FALSE;
    if (!checkDimensions)
    {
      std::cout << "\nError: number of dimensions is inconsistent.\n";
      std::cout << "Check boundary conditions in flowDescription.yaml.\n" << std::endl;
      exit(0);
    }
    BoundaryLocation location;
    for (PetscInt i=0; i<bcs.size(); i++)
    {
      location = stringToBoundaryLocation(bcs[i]["location"].as<std::string>());
      boundaries[location][U].type = stringToBoundaryType(bcs[i]["u"][0].as<std::string>());
      boundaries[location][U].value = bcs[i]["u"][1].as<PetscReal>();
      boundaries[location][V].type = stringToBoundaryType(bcs[i]["v"][0].as<std::string>());
      boundaries[location][V].value = bcs[i]["v"][1].as<PetscReal>();
      if (dim == 3)
      {
        boundaries[location][W].type = stringToBoundaryType(bcs[i]["w"][0].as<std::string>());
        boundaries[location][W].value = bcs[i]["w"][1].as<PetscReal>();
      }
    }
    
    // run some sanity checks on the input data
    PetscBool error = PETSC_FALSE;
    // if the boundary condition is periodic on one side, it should be periodic on the other side
    if (boundaries[XMINUS][U].type == PERIODIC && boundaries[XPLUS][U].type != PERIODIC)
      error = PETSC_TRUE;
    if (boundaries[XMINUS][U].type != PERIODIC && boundaries[XPLUS][U].type == PERIODIC)
      error = PETSC_TRUE;
    if (boundaries[YMINUS][U].type == PERIODIC && boundaries[YPLUS][U].type != PERIODIC)
      error = PETSC_TRUE;
    if (boundaries[YMINUS][U].type != PERIODIC && boundaries[YPLUS][U].type == PERIODIC)
      error = PETSC_TRUE;
    if (dim == 3)
    {
      if (boundaries[ZMINUS][U].type == PERIODIC && boundaries[ZPLUS][U].type != PERIODIC)
        error = PETSC_TRUE;
      if (boundaries[ZMINUS][U].type != PERIODIC && boundaries[ZPLUS][U].type == PERIODIC)
        error = PETSC_TRUE;
    }
    // if the boundary condition is periodic for one component, it should periodic for the others
    if (boundaries[XMINUS][U].type == PERIODIC && boundaries[XMINUS][V].type != PERIODIC)
      error = PETSC_TRUE;
    if (boundaries[XPLUS][U].type == PERIODIC && boundaries[XPLUS][V].type != PERIODIC)
      error = PETSC_TRUE;
    if (boundaries[YMINUS][U].type == PERIODIC && boundaries[YMINUS][V].type != PERIODIC)
      error = PETSC_TRUE;
    if (boundaries[YPLUS][U].type == PERIODIC && boundaries[YPLUS][V].type != PERIODIC)
      error = PETSC_TRUE;
    if (dim == 3)
    {
      if (boundaries[XMINUS][U].type == PERIODIC && boundaries[XMINUS][W].type != PERIODIC)
        error = PETSC_TRUE;
      if (boundaries[XPLUS][U].type == PERIODIC && boundaries[XPLUS][W].type != PERIODIC)
        error = PETSC_TRUE;
      if (boundaries[YMINUS][U].type == PERIODIC && boundaries[YMINUS][W].type != PERIODIC)
        error = PETSC_TRUE;
      if (boundaries[YPLUS][U].type == PERIODIC && boundaries[YPLUS][W].type != PERIODIC)
        error = PETSC_TRUE;
      if (boundaries[ZMINUS][U].type == PERIODIC && boundaries[ZMINUS][V].type != PERIODIC)
        error = PETSC_TRUE;
      if (boundaries[ZMINUS][U].type == PERIODIC && boundaries[ZMINUS][W].type != PERIODIC)
        error = PETSC_TRUE;
      if (boundaries[ZPLUS][U].type == PERIODIC && boundaries[ZPLUS][V].type != PERIODIC)
        error = PETSC_TRUE;
      if (boundaries[ZPLUS][U].type == PERIODIC && boundaries[ZPLUS][W].type != PERIODIC)
        error = PETSC_TRUE;
    }
    if (error)
    {
      std::cout << "\nERROR: Boundary conditions are inconsistent.\n";
      std::cout << "Check boundary conditions in flowDescription.yaml.\n" << std::endl;
      exit(0);
    }
  }
  MPI_Barrier(PETSC_COMM_WORLD);
  
  // broadcast flow description to all processes
  MPI_Bcast(&nu, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
  MPI_Bcast(initialVelocity, dim, MPIU_REAL, 0, PETSC_COMM_WORLD);
  MPI_Bcast(&initialCustomField, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
  MPI_Bcast(&perturbationAmplitude, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
  MPI_Bcast(&perturbationFrequency, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
  
  // create custom MPI type to broadcast boundary information
  MPI_Datatype boundaryInfoType, types[2];
  PetscInt blockcounts[2];
  MPI_Aint offsets[2];
  offsets[0] = offsetof(Boundary, type);
  types[0] = MPIU_INT;
  blockcounts[0] = 1;
  offsets[1] = offsetof(Boundary, value);
  types[1] = MPIU_REAL;
  blockcounts[1] = 1;
  MPI_Type_create_struct(2, blockcounts, offsets, types, &boundaryInfoType);
  MPI_Type_commit(&boundaryInfoType);
  MPI_Bcast(boundaries[XMINUS], dim, boundaryInfoType, 0, PETSC_COMM_WORLD);
  MPI_Bcast(boundaries[XPLUS], dim, boundaryInfoType, 0, PETSC_COMM_WORLD);
  MPI_Bcast(boundaries[YMINUS], dim, boundaryInfoType, 0, PETSC_COMM_WORLD);
  MPI_Bcast(boundaries[YPLUS], dim, boundaryInfoType, 0, PETSC_COMM_WORLD);
  if (dim == 3)
  {
    MPI_Bcast(boundaries[ZMINUS], dim, boundaryInfoType, 0, PETSC_COMM_WORLD);
    MPI_Bcast(boundaries[ZPLUS], dim, boundaryInfoType, 0, PETSC_COMM_WORLD);
  }
  MPI_Type_free(&boundaryInfoType);

  PetscPrintf(PETSC_COMM_WORLD, "done.\n");

} // initialize


/**
 * \brief Prints info about the initial and boundary conditions of the flow.
 */
template <PetscInt dim>
PetscErrorCode FlowDescription<dim>::printInfo()
{
  PetscErrorCode ierr;
  ierr = PetscPrintf(PETSC_COMM_WORLD, "\n---------------------------------------\n"); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Flow\n"); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "---------------------------------------\n"); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "dimensions: %d\n", dim); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "viscosity: %g\n", nu); CHKERRQ(ierr);
  if (initialCustomField)
  {
    ierr = PetscPrintf(PETSC_COMM_WORLD, "using an initial custom field\n"); CHKERRQ(ierr);
  }
  else
  {
    ierr = PetscPrintf(PETSC_COMM_WORLD, "initial velocity field:\n"); CHKERRQ(ierr);
    for (PetscInt i=0; i<dim; i++)
    {
      ierr = PetscPrintf(PETSC_COMM_WORLD, "\t%g\n", initialVelocity[i]); CHKERRQ(ierr);
    }
  }
  if (perturbationAmplitude != 0.0)
  {
    ierr = PetscPrintf(PETSC_COMM_WORLD, "using an initial perturbation:\n"); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "\tamplitude: %g\n", perturbationAmplitude); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "\tfrequency: %g\n", perturbationFrequency); CHKERRQ(ierr);
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD, "boundary conditions (component, type, value):\n"); CHKERRQ(ierr);
  for (PetscInt i=0; i<2*dim; i++)
  {
    ierr = PetscPrintf(PETSC_COMM_WORLD, "\t->location: %s\n", (stringFromBoundaryLocation(static_cast<BoundaryLocation>(i))).c_str()); CHKERRQ(ierr);
    for (PetscInt j=0; j<dim; j++)
    {
      ierr = PetscPrintf(PETSC_COMM_WORLD, "\t\t%d \t %s \t %g\n", 
                                           j, 
                                           (stringFromBoundaryType(boundaries[i][j].type).c_str()),
                                           boundaries[i][j].value); CHKERRQ(ierr);
    }
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD, "---------------------------------------\n"); CHKERRQ(ierr);

  return 0;
} // printInfo


// dimension spacialization
template class FlowDescription<2>;
template class FlowDescription<3>;