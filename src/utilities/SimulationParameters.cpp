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
      exit(0);
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
      exit(0);
      break;
  }

  const YAML::Node &solvers = node["linearSolvers"];
  for (unsigned int i=0; i<solvers.size(); i++)
  {
    if (solvers[i]["system"].as<std::string>() == "velocity")
    {
      velocitySolver.method = stringToIterativeMethod(solvers[i]["solver"].as<std::string>("CG"));
      velocitySolver.preconditioner = stringToPreconditionerType(solvers[i]["preconditioner"].as<std::string>("DIAGONAL"));
      velocitySolver.relativeTolerance = solvers[i]["relativeTolerance"].as<PetscReal>(1.0E-05);
      velocitySolver.maxIterations = solvers[i]["maxIterations"].as<PetscInt>(10000);
    }
    else if (solvers[i]["system"].as<std::string>() == "Poisson" || solvers[i]["system"].as<std::string>() == "poisson")
    {
      poissonSolver.method = stringToIterativeMethod(solvers[i]["solver"].as<std::string>("CG"));
      poissonSolver.preconditioner = stringToPreconditionerType(solvers[i]["preconditioner"].as<std::string>("DIAGONAL"));
      poissonSolver.relativeTolerance = solvers[i]["relativeTolerance"].as<PetscReal>(1.0E-05);
      poissonSolver.maxIterations = solvers[i]["maxIterations"].as<PetscInt>(10000);
    }
    else
    {
      std::cout << "\nERROR: " << solvers[i]["system"].as<std::string>() << " - unknown solver name.\n";
      std::cout << "solver name implemented:\n";
      std::cout << "\tvelocity (solver for the intermediate velocity)\n";
      std::cout << "\tpoisson (poisson solver, modified or not)\n" << std::endl;
      exit(0);
    }
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
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Time-stepping\n"); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "---------------------------------------\n"); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "formulation: %s\n", stringFromIBMethod(ibm).c_str()); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "convection: %s\n", stringFromTimeScheme(convection.scheme).c_str()); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "diffusion: %s\n", stringFromTimeScheme(diffusion.scheme).c_str()); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "time-increment: %g\n", dt); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "starting time-step: %d\n", startStep); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "number of time-steps: %d\n", nt); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "saving-interval: %d\n", nsave); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "---------------------------------------\n"); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD, "\n---------------------------------------\n"); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Linear system for intermediate velocity\n"); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "---------------------------------------\n"); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "method: %s\n", stringFromIterativeMethod(velocitySolver.method).c_str()); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "preconditioner: %s\n", stringFromPreconditionerType(velocitySolver.preconditioner).c_str()); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "relative-tolerance: %g\n", velocitySolver.relativeTolerance); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "maximum number of iterations: %d\n", velocitySolver.maxIterations); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "---------------------------------------\n"); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD, "\n---------------------------------------\n"); CHKERRQ(ierr);
  switch (ibm)
  {
    case NAVIER_STOKES:
      ierr = PetscPrintf(PETSC_COMM_WORLD, "Linear system for pressure\n"); CHKERRQ(ierr);
      break;
    case TAIRA_COLONIUS:
      ierr = PetscPrintf(PETSC_COMM_WORLD, "Linear system for pressure-force\n"); CHKERRQ(ierr);
      break;
    default:
      break;
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD, "---------------------------------------\n"); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "method: %s\n", stringFromIterativeMethod(poissonSolver.method).c_str()); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "preconditioner: %s\n", stringFromPreconditionerType(poissonSolver.preconditioner).c_str()); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "relative-tolerance: %g\n", poissonSolver.relativeTolerance); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "maximum number of iterations: %d\n", poissonSolver.maxIterations); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "---------------------------------------\n"); CHKERRQ(ierr);

  return 0;
} // printInfo