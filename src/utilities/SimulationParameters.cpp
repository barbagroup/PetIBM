/***************************************************************************//**
 * \file SimulationParameters.cpp
 * \author Anush Krishnan (anush@bu.edu), Olivier Mesnard (mesnardo@gwu.edu)
 * \brief Implementation of the methods of the class \c SimulationParameters.
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
 * \brief Constructor -- Parses simulationParameters.yaml.
 */
SimulationParameters::SimulationParameters(std::string directory)
{
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
 * \param filePath path of the simulation parameters file
 */
void SimulationParameters::initialize(std::string filePath)
{
  PetscPrintf(PETSC_COMM_WORLD, "\nParsing file %s... ", filePath.c_str());

  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank); // get rank of current process
  
  if (rank == 0) // read the input file only on process 0
  {
    YAML::Node nodes(YAML::LoadFile(filePath));
    const YAML::Node &node = nodes[0];

    dt = node["dt"].as<PetscReal>();
    startStep = node["startStep"].as<PetscInt>(0);
    nt = node["nt"].as<PetscInt>();
    nsave = node["nsave"].as<PetscInt>(nt);

    ibmScheme = stringToIBMScheme(node["ibmScheme"].as<std::string>("NAVIER_STOKES"));

    convection.scheme = stringToTimeScheme(node["convection"].as<std::string>("EULER_EXPLICIT"));
    diffusion.scheme = stringToTimeScheme(node["diffusion"].as<std::string>("EULER_IMPLICIT"));
    // set time-stepping coefficients for convective terms
    switch (convection.scheme)
    {
      case EULER_EXPLICIT:
        convection.coefficients.push_back(0.0); // n+1 coefficient
        convection.coefficients.push_back(1.0); // n coefficient
        break;
      case ADAMS_BASHFORTH_2:
        convection.coefficients.push_back(1.5);  // n+1 coefficient
        convection.coefficients.push_back(-0.5); // n coefficient
        break;
      default:
        std::cout << "\nERROR: unknown numerical scheme for convective terms.\n";
        std::cout << "Numerical scheme implemented:\n";
        std::cout << "\tEULER_EXPLICIT\n";
        std::cout << "\tADAMS_BASHFORTH_2\n" << std::endl;
        exit(0);
        break;
    }
    // set time-stepping coefficients for diffusive terms
    switch (diffusion.scheme)
    {
      case EULER_EXPLICIT:
        diffusion.coefficients.push_back(0.0); // n+1 coefficient
        diffusion.coefficients.push_back(1.0); // n coefficient
        break;
      case EULER_IMPLICIT:
        diffusion.coefficients.push_back(1.0); // n+1 coefficient
        break;
      case CRANK_NICOLSON:
        diffusion.coefficients.push_back(0.5); // n+1 coefficient
        diffusion.coefficients.push_back(0.5); // n coefficient
        break;
      default:
        std::cout << "\nERROR: unknown numerical scheme for diffusive terms.\n";
        std::cout << "Numerical scheme implemented:\n";
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
  }
  MPI_Barrier(PETSC_COMM_WORLD);
  
  // broadcast parameters to all processes
  MPI_Bcast(&dt, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
  MPI_Bcast(&nt, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
  MPI_Bcast(&nsave, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
  MPI_Bcast(&startStep, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
  
  MPI_Bcast(&ibmScheme, 1, MPI_CHAR, 0, PETSC_COMM_WORLD);

  // create custom MPI type to broadcast time-stepping schemes
  MPI_Datatype timeIntegrationInfoType;
  MPI_Datatype types1[2] = {MPIU_INT, MPIU_REAL};
  PetscInt blockcounts1[2];
  MPI_Aint offsets1[2];
  offsets1[0] = offsetof(TimeIntegration, scheme);
  offsets1[1] = offsetof(TimeIntegration, coefficients);
  blockcounts1[0] = 1;
  // broadcast convection scheme
  blockcounts1[1] = convection.coefficients.size();
  MPI_Type_create_struct(2, blockcounts1, offsets1, types1, &timeIntegrationInfoType);
  MPI_Type_commit(&timeIntegrationInfoType);
  MPI_Bcast(&convection, 1, timeIntegrationInfoType, 0, PETSC_COMM_WORLD);
  // broadcast diffusion scheme
  blockcounts1[1] = diffusion.coefficients.size();
  MPI_Type_create_struct(2, blockcounts1, offsets1, types1, &timeIntegrationInfoType);
  MPI_Type_commit(&timeIntegrationInfoType);
  MPI_Bcast(&diffusion, 1, timeIntegrationInfoType, 0, PETSC_COMM_WORLD);
  MPI_Type_free(&timeIntegrationInfoType);

  // create custom MPI type to broadcast solver information
  MPI_Datatype solverInfoType;
  MPI_Datatype types2[4] = {MPIU_INT, MPIU_INT, MPIU_REAL, MPI_INT};
  PetscInt blockcounts2[4] = {1, 1, 1, 1};
  MPI_Aint offsets2[4];
  offsets2[0] = offsetof(Solver, method);
  offsets2[1] = offsetof(Solver, preconditioner);
  offsets2[2] = offsetof(Solver, relativeTolerance);
  offsets2[3] = offsetof(Solver, maxIterations);
  MPI_Type_create_struct(4, blockcounts2, offsets2, types2, &solverInfoType);
  MPI_Type_commit(&solverInfoType);
  MPI_Bcast(&velocitySolver, 1, solverInfoType, 0, PETSC_COMM_WORLD);
  MPI_Bcast(&poissonSolver, 1, solverInfoType, 0, PETSC_COMM_WORLD);
  MPI_Type_free(&solverInfoType);

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
  ierr = PetscPrintf(PETSC_COMM_WORLD, "formulation: %s\n", stringFromIBMScheme(ibmScheme).c_str()); CHKERRQ(ierr);
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
  switch (ibmScheme)
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