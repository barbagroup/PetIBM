/***************************************************************************//**
 * \file SimulationParameters.cpp
 * \author Anush Krishnan (anush@bu.edu), Olivier Mesnard (mesnardo@gwu.edu)
 * \brief Implementation of the methods of the class `SimulationParameters`.
 */


#include "SimulationParameters.h"

#include <fstream>

#include "yaml-cpp/yaml.h"


/**
 * \brief Constructor.
 */
SimulationParameters::SimulationParameters()
{
} // SimulationParameters


/**
 * \brief Constructor -- Parses the input file `simulationParameters.yaml`.
 *
 * \param dir Directory of the simulation
 */
SimulationParameters::SimulationParameters(std::string dir)
{
  directory = dir;
  initialize(directory + "/simulationParameters.yaml");
} // SimulationParameters


/**
 * \brief Destructor.
 */
SimulationParameters::~SimulationParameters()
{
} // ~SimulationParameters


/**
 * \brief Parses file containing the simulation parameters.
 *
 * The file is parsed using YAML format.
 *
 * \param filePath Path of the simulation parameters file
 */
void SimulationParameters::initialize(std::string filePath)
{
  PetscPrintf(PETSC_COMM_WORLD, "\nParsing file %s... ", filePath.c_str());
  
  YAML::Node nodes(YAML::LoadFile(filePath));
  const YAML::Node &node = nodes[0];

  dt = node["dt"].as<PetscReal>();
  startStep = node["startStep"].as<PetscInt>(0);
  nt = node["nt"].as<PetscInt>();
  nsave = node["nsave"].as<PetscInt>(nt);
  
  vSolveType = stringToExecuteType(node["vSolveType"].as<std::string>("CPU"));
  pSolveType = stringToExecuteType(node["pSolveType"].as<std::string>("CPU"));

  outputFormat = node["outputFormat"].as<std::string>("binary");
#ifndef PETSC_HAVE_HDF5
  if (outputFormat == "hdf5")
  {
    PetscPrintf(PETSC_COMM_WORLD,
                "\nERROR: PETSc has not been built with HDF5 available; "
                "you cannot use `outputFormat: hdf5`\n");
    MPI_Barrier(PETSC_COMM_WORLD);
    exit(1);
  }
#endif
  outputFlux = (node["outputFlux"].as<bool>(true)) ? PETSC_TRUE : PETSC_FALSE;
  outputVelocity = (node["outputVelocity"].as<bool>(false)) ? PETSC_TRUE : PETSC_FALSE;

  ibm = stringToIBMethod(node["ibm"].as<std::string>("NONE"));

  convection.scheme = stringToTimeScheme(node["convection"].as<std::string>("EULER_EXPLICIT"));
  diffusion.scheme = stringToTimeScheme(node["diffusion"].as<std::string>("EULER_IMPLICIT"));
  // set time-stepping coefficients for convective terms
  switch (convection.scheme)
  {
    case NONE:
      convection.coefficients.push_back(0.0); // n+1 coefficient
      convection.coefficients.push_back(0.0); // n coefficient
      convection.coefficients.push_back(0.0); // n-1 coefficient
      break;
    case EULER_EXPLICIT:
      convection.coefficients.push_back(0.0); // n+1 coefficient
      convection.coefficients.push_back(1.0); // n coefficient
      convection.coefficients.push_back(0.0); // n-1 coefficient
      break;
    case ADAMS_BASHFORTH_2:
      convection.coefficients.push_back(0.0);  // n+1 coefficient
      convection.coefficients.push_back(1.5);  // n coefficient
      convection.coefficients.push_back(-0.5); // n-1 coefficient
      break;
    default:
      std::cout << "\nERROR: unknown numerical scheme for convective terms.\n";
      std::cout << "Numerical scheme implemented:\n";
      std::cout << "\tNONE\n";
      std::cout << "\tEULER_EXPLICIT\n";
      std::cout << "\tADAMS_BASHFORTH_2\n" << std::endl;
      exit(1);
      break;
  }
  // set time-stepping coefficients for diffusive terms
  switch (diffusion.scheme)
  {
    case NONE:
      diffusion.coefficients.push_back(0.0); // n+1 coefficient
      diffusion.coefficients.push_back(0.0); // n coefficient
      break;
    case EULER_EXPLICIT:
      diffusion.coefficients.push_back(0.0); // n+1 coefficient
      diffusion.coefficients.push_back(1.0); // n coefficient
      break;
    case EULER_IMPLICIT:
      diffusion.coefficients.push_back(1.0); // n+1 coefficient
      diffusion.coefficients.push_back(0.0); // n coefficient
      break;
    case CRANK_NICOLSON:
      diffusion.coefficients.push_back(0.5); // n+1 coefficient
      diffusion.coefficients.push_back(0.5); // n coefficient
      break;
    default:
      std::cout << "\nERROR: unknown numerical scheme for diffusive terms.\n";
      std::cout << "Numerical scheme implemented:\n";
      std::cout << "\tNONE\n";
      std::cout << "\tEULER_EXPLICIT\n";
      std::cout << "\tEULER_IMPLICIT\n";
      std::cout << "\tCRANK_NICOLSON\n" << std::endl;
      exit(1);
      break;
  }

  PetscPrintf(PETSC_COMM_WORLD, "done.\n");
  
} // initialize


/**
 * \brief Prints info about the initial and boundary conditions of the flow.
 */
PetscErrorCode SimulationParameters::printInfo()
{
  PetscErrorCode ierr;
  ierr = PetscPrintf(PETSC_COMM_WORLD, "\n---------------------------------------\n"); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Simulation parameters\n"); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "---------------------------------------\n"); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "formulation: %s\n", stringFromIBMethod(ibm).c_str()); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "convection: %s\n", stringFromTimeScheme(convection.scheme).c_str()); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "diffusion: %s\n", stringFromTimeScheme(diffusion.scheme).c_str()); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "time-increment: %g\n", dt); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "starting time-step: %d\n", startStep); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "number of time-steps: %d\n", nt); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "saving-interval: %d\n", nsave); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "velocity solver type: %s\n", stringFromExecuteType(vSolveType).c_str()); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Poisson solver type: %s\n", stringFromExecuteType(pSolveType).c_str()); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "output format: %s\n", outputFormat.c_str()); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "output flux: %D\n", outputFlux); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "output velocity: %D\n", outputVelocity); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "---------------------------------------\n"); CHKERRQ(ierr);

  return 0;
} // printInfo