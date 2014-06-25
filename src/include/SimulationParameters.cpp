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

/***************************************************************************//**
* \param fileName Input file path
*
* This is the constructor for the class SimulationParameters. A case folder with
* the input files is supplied to the flow solver, and this function reads the
* file \c simulationParameters.yaml in the folder. The parameters are listed in
* the file using the YAML format.
*/
SimulationParameters::SimulationParameters(std::string fileName)
{
	PetscInt    rank;
	
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank); // get rank of current process
	
	if(rank==0) // read the input file only on process 0
	{
		std::ifstream inputFile(fileName.c_str());
		YAML::Parser  parser(inputFile);
		YAML::Node    document;
		std::string   solverString;
		std::string   convSchString, diffSchString;
		std::string   restartString = "no";
		
		parser.GetNextDocument(document);

		document[0]["dt"] >> dt;       // read the time step size
		document[0]["nt"] >> nt;       // read the number of time steps
		document[0]["nsave"] >> nsave; // read the save interval
		
		restart = PETSC_FALSE;         // assume no restart by default
		startStep = 0;                 // the simulation starts from time step 0

		// read restartString from the input file
		// this input is optional, and hence is inside a try block
		// the default value of restartString (set above) is "no"
		try
		{
			document[0]["restart"] >> restartString;
		}
		catch(...)
		{
		}
		// set restart depending on the value of restartString
		if(restartString=="yes" || restartString=="true") restart = PETSC_TRUE;
		if(restart)
		{
			// read the starting time step for the restarted simulation
			try
			{
				document[0]["startStep"] >> startStep;
			}
			catch(...)
			{
			}
		}
		document[0]["ibmScheme"] >> solverString;      // read the type of flow solver
		document[0]["timeScheme"][0] >> convSchString; // read the time-stepping scheme for convection
		document[0]["timeScheme"][1] >> diffSchString; // read the time-stepping scheme for diffusion
		
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
		inputFile.close();
	}
	MPI_Barrier(PETSC_COMM_WORLD);
	
	// broadcast parameters to all processes
	MPI_Bcast(&dt, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
	MPI_Bcast(&nt, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
	MPI_Bcast(&nsave, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
	MPI_Bcast(&restart, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
	MPI_Bcast(&startStep, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
	MPI_Bcast(&gamma, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
	MPI_Bcast(&zeta, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
	MPI_Bcast(&alphaExplicit, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
	MPI_Bcast(&alphaImplicit, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
	MPI_Bcast(&solverType, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
	MPI_Bcast(&convectionScheme, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
	MPI_Bcast(&diffusionScheme, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
}
