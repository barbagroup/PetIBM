/***************************************************************************//**
 * \file SimulationParameters.cpp
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the methods of the class \c SimulationParameters.
 */


#include "SimulationParameters.h"

#include <fstream>

#include "yaml-cpp/yaml.h"


/**
 * \brief Converts \c std::string to \c TimeSteppingScheme.
 */
TimeSteppingScheme timeSchemeFromString(std::string &s)
{
  if (s == "EULER_EXPLICIT")
    return EULER_EXPLICIT;
  if (s == "EULER_IMPLICIT")
    return EULER_IMPLICIT;
  if (s == "ADAMS_BASHFORTH_2")
    return ADAMS_BASHFORTH_2;
  if (s == "CRANK_NICOLSON")
    return CRANK_NICOLSON;
  return EULER_EXPLICIT;
}

/**
 * \brief Converts \c std::string to \c SolverType.
 */
SolverType solverTypeFromString(std::string &s)
{
  if (s == "NAVIER_STOKES")
    return NAVIER_STOKES;
  if (s == "TAIRA_COLONIUS")
    return TAIRA_COLONIUS;
  return NAVIER_STOKES;
}

SimulationParameters::SimulationParameters()
{
}

/**
 * \brief Constructor -- Parses simulationParameters.yaml
 */
SimulationParameters::SimulationParameters(std::string fileName)
{
  initialize(fileName);
}

/**
 * \brief Parses file containing the simulation parameters.
 *
 * The file is parsed using YAML format.
 *
 * \param fileName path of the simulation parameters file
 */
void SimulationParameters::initialize(std::string fileName)
{
  PetscInt    rank;
  
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank); // get rank of current process
  
  if (rank == 0) // read the input file only on process 0
  {
    YAML::Node nodes(YAML::LoadFile(fileName));

    dt = nodes[0]["dt"].as<PetscReal>();
    startStep = nodes[0]["startStep"].as<PetscInt>(0);
    restart = (startStep > 0) ? PETSC_TRUE : PETSC_FALSE;
    nt = nodes[0]["nt"].as<PetscInt>();
    nsave = nodes[0]["nsave"].as<PetscInt>(nt);

    solverType = solverTypeFromString(nodes[0]["ibmScheme"].as<std::string>("NAVIER_STOKES"));
    convectionScheme = timeSchemeFromString(nodes[0]["timeScheme"][0].as<std::string>("EULER_EXPLICIT"));
    diffusionScheme  = timeSchemeFromString(nodes[0]["timeScheme"][1].as<std::string>("EULER_IMPLICIT"));

    // set the time-stepping coefficients for the different schemes
    switch (convectionScheme)
    {
      case EULER_EXPLICIT:
        gamma = 1.0;
        zeta = 0.0;
        break;
      case ADAMS_BASHFORTH_2:
        gamma = 1.5;
        zeta = -0.5;
        break;
      default:
        gamma = 1.0;
        zeta = 0.0;
        break;
    }
    switch (diffusionScheme)
    {
      case EULER_EXPLICIT:
        alphaExplicit = 1.0;
        alphaImplicit = 0.0;
        break;
      case EULER_IMPLICIT:
        alphaExplicit = 0.0;
        alphaImplicit = 1.0;
        break;
      case CRANK_NICOLSON:
        alphaExplicit = 0.5;
        alphaImplicit = 0.5;
        break;
      default:
        alphaExplicit = 1.0;
        alphaImplicit = 0.0;
        break;
    }

    const YAML::Node &systems = nodes[0]["linearSystems"];
    std::string name, solver, preconditioner;
    for (unsigned int i=0; i<systems.size(); i++)
    {
      name = systems[i]["system"].as<std::string>();
      solver = systems[i]["solver"].as<std::string>("CG");
      preconditioner = systems[i]["preconditioner"].as<std::string>("DIAGONAL");

      // set the simulation parameters
      if (name == "velocity")
      {
        velocitySolveTolerance = systems[i]["tolerance"].as<PetscReal>(1.0E-05);
        velocitySolveMaxIts = systems[i]["maxIterations"].as<PetscInt>(10000);
      }
      else if (name == "Poisson")
      {
        PoissonSolveTolerance = systems[i]["tolerance"].as<PetscReal>(1.0E-05);
        PoissonSolveMaxIts = systems[i]["maxIterations"].as<PetscInt>(10000);
      }
    }
  }
  MPI_Barrier(PETSC_COMM_WORLD);
  
  // broadcast parameters to all processes
  MPI_Bcast(&dt, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
  MPI_Bcast(&nt, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
  MPI_Bcast(&nsave, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
  MPI_Bcast(&startStep, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
  MPI_Bcast(&restart, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
  
  MPI_Bcast(&gamma, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
  MPI_Bcast(&zeta, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
  MPI_Bcast(&alphaExplicit, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
  MPI_Bcast(&alphaImplicit, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
  
  MPI_Bcast(&solverType, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
  
  MPI_Bcast(&convectionScheme, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
  MPI_Bcast(&diffusionScheme, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
  
  MPI_Bcast(&velocitySolveTolerance, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
  MPI_Bcast(&PoissonSolveTolerance, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
  MPI_Bcast(&velocitySolveMaxIts, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
  MPI_Bcast(&PoissonSolveMaxIts, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
}