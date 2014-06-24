#include "SimulationParameters.h"
#include "yaml-cpp/yaml.h"
#include <fstream>

TimeSteppingScheme timeSchemeFromString(std::string &s)
{
	if (s == "EULER_EXPLICIT") return EULER_EXPLICIT;
	if (s == "EULER_IMPLICIT") return EULER_IMPLICIT;
	if (s == "ADAMS_BASHFORTH_2") return ADAMS_BASHFORTH_2;
	if (s == "RUNGE_KUTTA_3") return RUNGE_KUTTA_3;
	if (s == "CRANK_NICOLSON") return CRANK_NICOLSON;
	return EULER_EXPLICIT;
}

SolverType solverTypeFromString(std::string &s)
{
	if (s == "NAVIER_STOKES") return NAVIER_STOKES;
	if (s == "SAIKI_BIRINGEN") return SAIKI_BIRINGEN;
	if (s == "FADLUN_ET_AL") return FADLUN_ET_AL;
	if (s == "TAIRA_COLONIUS") return TAIRA_COLONIUS;
	return NAVIER_STOKES;
}

SimulationParameters::SimulationParameters(std::string fileName)
{
	PetscInt    rank;
	
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	if(rank==0)
	{
		std::ifstream file(fileName.c_str());
		YAML::Parser  parser(file);
		YAML::Node    doc;
		std::string   solver;
		std::string   convSch, diffSch;
		std::string   restartOption = "no";
		
		parser.GetNextDocument(doc);
				
		doc[0]["dt"] >> dt;
		doc[0]["nt"] >> nt;
		doc[0]["nsave"] >> nsave;
		restart = PETSC_FALSE;
		startStep = 0;
		try
		{
			doc[0]["restart"] >> restartOption;
		}
		catch(...)
		{
		}
		if(restartOption == "yes" || restartOption == "true") restart = PETSC_TRUE;
		if(restart)
		{
			try
			{
				doc[0]["startStep"] >> startStep;
			}
			catch(...)
			{
			}
		}
		doc[0]["ibmScheme"] >> solver;
		doc[0]["timeScheme"][0] >> convSch;
		doc[0]["timeScheme"][1] >> diffSch;
		convectionScheme = timeSchemeFromString(convSch);
		diffusionScheme  = timeSchemeFromString(diffSch);
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
		solverType = solverTypeFromString(solver);
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
