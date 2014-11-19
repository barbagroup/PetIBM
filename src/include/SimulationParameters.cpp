/***************************************************************************//**
* \file
* \brief Source file to define member functions of SimulationParameters
*/

#include "SimulationParameters.h"
#include "yaml-cpp/yaml.h"
#include <fstream>

/***************************************************************************//**
* \brief Converts \c std::string to \c TimeSteppingScheme
*/
TimeSteppingScheme timeSchemeFromString(std::string &s)
{
	if (s == "EULER_EXPLICIT") return EULER_EXPLICIT;
	if (s == "EULER_IMPLICIT") return EULER_IMPLICIT;
	if (s == "ADAMS_BASHFORTH_2") return ADAMS_BASHFORTH_2;
	if (s == "CRANK_NICOLSON") return CRANK_NICOLSON;
	return EULER_EXPLICIT;
}

/***************************************************************************//**
* \brief Converts \c std::string to \c SolverType.
*/
SolverType solverTypeFromString(std::string &s)
{
	if (s == "NAVIER_STOKES") return NAVIER_STOKES;
	if (s == "TAIRA_COLONIUS") return TAIRA_COLONIUS;
	return NAVIER_STOKES;
}

SimulationParameters::SimulationParameters()
{
}

SimulationParameters::SimulationParameters(std::string fileName)
{
	initialize(fileName);
}

/***************************************************************************//**
* \param fileName Input file path
*
* This is the constructor for the class SimulationParameters. A case folder with
* the input files is supplied to the flow solver, and this function reads the
* file \c simulationParameters.yaml in the folder. The parameters are listed in
* the file using the YAML format.
*/
void SimulationParameters::initialize(std::string fileName)
{
	PetscInt    rank;
	
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank); // get rank of current process
	
	if(rank==0) // read the input file only on process 0
	{
		std::ifstream inputFile(fileName.c_str());
		YAML::Parser  parser(inputFile);
		YAML::Node    document;
		std::string   solverString  = "NAVIER_STOKES";
		std::string   convSchString = "EULER_EXPLICIT",
		              diffSchString = "EULER_IMPLICIT";
		std::string   restartString = "no";
		
		parser.GetNextDocument(document);

		document[0]["dt"] >> dt;       // read the time step size
		document[0]["nt"] >> nt;       // read the number of time steps
		document[0]["nsave"] >> nsave; // read the save interval
		
		startStep = 0;                 // the simulation starts from time step 0
		restart = PETSC_FALSE;         // assume no restart by default

		// read the starting time step for the restarted simulation
		// this input is optional, and hence is inside a try block
		// the default value of startStep (set above) is 0
		try
		{
			document[0]["startStep"] >> startStep;
		}
		catch(...)
		{
		}

		// set restart to true of startStep is different from 0
		if(startStep > 0) restart = PETSC_TRUE;

		document[0]["ibmScheme"] >> solverString;      // read the type of flow solver

		try
		{
			document[0]["timeScheme"][0] >> convSchString; // read the time-stepping scheme for convection
			document[0]["timeScheme"][1] >> diffSchString; // read the time-stepping scheme for diffusion
		}
		catch(...)
		{
		}
		
		// set the solver type and time-stepping schemes
		// from the strings read from the input file
		solverType       = solverTypeFromString(solverString);
		convectionScheme = timeSchemeFromString(convSchString);
		diffusionScheme  = timeSchemeFromString(diffSchString);
		
		// set the time-stepping coefficients for the different schemes
		switch(convectionScheme)
		{
			case EULER_EXPLICIT:
				gamma = 1.0;
				zeta  = 0.0;
				break;
			case ADAMS_BASHFORTH_2:
				gamma = 1.5;
				zeta  = -0.5;
				break;
			default:
				gamma = 1.0;
				zeta  = 0.0;
				break;
		}
		switch(diffusionScheme)
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

		std::string system = "",
		            linearSolver = "",
		            preconditioner = "";
		PetscReal   tol = 1e-5;
		PetscInt    maxIter = 10000;

		velocitySolveTolerance = 1e-5;
		velocitySolveMaxIts    = 10000;
		PoissonSolveTolerance  = 1e-5;
		PoissonSolveMaxIts     = 20000;
		const YAML::Node &solvers = document[0]["linearSolvers"];
		for (unsigned int i=0; i<solvers.size(); i++)
		{
			// read linear solver options
			solvers[i]["system"] >> system;
			solvers[i]["solver"] >> linearSolver;
			try
			{
				solvers[i]["preconditioner"] >> preconditioner;
			}
			catch(...)
			{
			}
			try
			{
				solvers[i]["tolerance"] >> tol;
			}
			catch(...)
			{
			}
			try
			{
				solvers[i]["maxIterations"] >> maxIter;
			}
			catch(...)
			{
			}
			// set the simulation parameters
			if(system=="velocity")
			{
				velocitySolveTolerance = tol;
				velocitySolveMaxIts    = maxIter;
			}
			if(system=="Poisson")
			{
				PoissonSolveTolerance = tol;
				PoissonSolveMaxIts    = maxIter;
			}
		}
		inputFile.close();
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
