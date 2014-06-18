#include "FlowDescription.h"
#include "yaml-cpp/yaml.h"
#include <iostream>
#include <fstream>

Boundary boundaryFromString(std::string &s)
{
	if (s == "xMinus") return XMINUS;
	if (s == "xPlus") return XPLUS;
	if (s == "yMinus") return YMINUS;
	if (s == "yPlus") return YPLUS;
	if (s == "zMinus") return ZMINUS;
	if (s == "zPlus") return ZPLUS;
	
	std::cout << "ERROR: Invalid boundary location!\n";
	exit(0);
}

BCType bcTypeFromString(std::string &s)
{
	if (s == "DIRICHLET") return DIRICHLET;
	if (s == "NEUMANN") return NEUMANN;
	if (s == "CONVECTIVE") return CONVECTIVE;
	if (s == "PERIODIC") return PERIODIC;
	
	std::cout << "ERROR: Invalid boundary condition type!\n";
	exit(0);
}

FlowDescription::FlowDescription(std::string fileName)
{
	PetscInt    rank;
	
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	if(rank==0)
	{
		std::ifstream file(fileName.c_str());
		YAML::Parser  parser(file);
		YAML::Node    doc;
		std::string   loc, typ;
		Boundary      location;
		
		parser.GetNextDocument(doc);
		
		doc[0]["dimensions"] >> dimensions;
		doc[0]["nu"] >> nu;
		
		doc[0]["initialVelocity"][0] >> initialVelocity[0];
		doc[0]["initialVelocity"][1] >> initialVelocity[1];
		if(dimensions==3)
			doc[0]["initialVelocity"][2] >> initialVelocity[2];
		else
			initialVelocity[2] = 0.0;
		
		try
		{
			doc[0]["initialPerturbation"][0] >> initialPerturbation[0];
			doc[0]["initialPerturbation"][1] >> initialPerturbation[1];
			if(dimensions==3)
				doc[0]["initialPerturbation"][2] >> initialPerturbation[2];
			else
				initialVelocity[2] = 0.0;
		}
		catch(...)
		{
		}
		
		const YAML::Node &BCs = doc[0]["boundaryConditions"];
		
		for (size_t i=0; i<BCs.size(); i++)
		{
			BCs[i]["location"] >> loc;
			location = boundaryFromString(loc);
			
			BCs[i]["u"][0] >> typ;
			bc[0][location].type = bcTypeFromString(typ);
			BCs[i]["u"][1] >> bc[0][location].value;
			
			BCs[i]["v"][0] >> typ;
			bc[1][location].type = bcTypeFromString(typ);
			BCs[i]["v"][1] >> bc[1][location].value;
			
			if(dimensions>2)
			{
				BCs[i]["w"][0] >> typ;
				bc[2][location].type = bcTypeFromString(typ);
				BCs[i]["w"][1] >> bc[2][location].value;
			}
			else
			{
				bc[2][location].type  = PERIODIC;
				bc[2][location].value = 0.0;
			}
		}
		
		PetscBool flag = PETSC_TRUE;
		// check if opposite faces are both periodic
		if(bc[0][XMINUS].type==PERIODIC && bc[0][XPLUS].type!=PERIODIC) flag = PETSC_FALSE;
		if(bc[0][XMINUS].type!=PERIODIC && bc[0][XPLUS].type==PERIODIC) flag = PETSC_FALSE;
		if(bc[0][YMINUS].type==PERIODIC && bc[0][YPLUS].type!=PERIODIC) flag = PETSC_FALSE;
		if(bc[0][YMINUS].type!=PERIODIC && bc[0][YPLUS].type==PERIODIC) flag = PETSC_FALSE;
		if(dimensions==3)
		{
			if(bc[0][ZMINUS].type==PERIODIC && bc[0][ZPLUS].type!=PERIODIC) flag = PETSC_FALSE;
			if(bc[0][ZMINUS].type!=PERIODIC && bc[0][ZPLUS].type==PERIODIC) flag = PETSC_FALSE;
		}
		// check if all components of velocity on the same face are periodic
		if(bc[0][XMINUS].type==PERIODIC && bc[1][XMINUS].type!=PERIODIC) flag = PETSC_FALSE;
		if(bc[0][XPLUS].type==PERIODIC && bc[1][XPLUS].type!=PERIODIC) flag = PETSC_FALSE;
		if(bc[0][YMINUS].type==PERIODIC && bc[1][YMINUS].type!=PERIODIC) flag = PETSC_FALSE;
		if(bc[0][YPLUS].type==PERIODIC && bc[1][YPLUS].type!=PERIODIC) flag = PETSC_FALSE;
		if(dimensions==3)
		{
			// check for w on X and Y faces
			if(bc[0][XMINUS].type==PERIODIC && bc[2][XMINUS].type!=PERIODIC) flag = PETSC_FALSE;
			if(bc[0][XPLUS].type==PERIODIC && bc[2][XPLUS].type!=PERIODIC) flag = PETSC_FALSE;
			if(bc[0][YMINUS].type==PERIODIC && bc[2][YMINUS].type!=PERIODIC) flag = PETSC_FALSE;
			if(bc[0][YPLUS].type==PERIODIC && bc[2][YPLUS].type!=PERIODIC) flag = PETSC_FALSE;
			
			// check for v and w on Z faces
			if(bc[0][ZMINUS].type==PERIODIC && bc[1][ZMINUS].type!=PERIODIC) flag = PETSC_FALSE;
			if(bc[0][ZMINUS].type==PERIODIC && bc[2][ZMINUS].type!=PERIODIC) flag = PETSC_FALSE;
			if(bc[0][ZPLUS].type==PERIODIC && bc[1][ZPLUS].type!=PERIODIC) flag = PETSC_FALSE;
			if(bc[0][ZPLUS].type==PERIODIC && bc[2][ZPLUS].type!=PERIODIC) flag = PETSC_FALSE;
		}
		if(!flag)
		{
			std::cout << "ERROR: Check if boundary conditions are consistent." << std::endl;
			exit(0);
		}
	}
	MPI_Barrier(PETSC_COMM_WORLD);
	
	// broadcast flow description to all processes
	MPI_Bcast(&dimensions, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
	MPI_Bcast(&nu, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
	MPI_Bcast(initialVelocity, 3, MPIU_REAL, 0, PETSC_COMM_WORLD);
	MPI_Bcast(initialPerturbation, 3, MPIU_REAL, 0, PETSC_COMM_WORLD);
	
	// create custom MPI type to broadcast BC information
	MPI_Datatype bcInfoType, types[2];
	PetscInt     blockcounts[2];
	MPI_Aint     offsets[2];
	offsets[0]     = offsetof(BoundaryCondition, type);
	types[0]       = MPIU_INT;
	blockcounts[0] = 1;
	offsets[1]     = offsetof(BoundaryCondition, value);
	types[1]       = MPIU_REAL;
	blockcounts[1] = 1;
	MPI_Type_create_struct(2, blockcounts, offsets, types, &bcInfoType);
	MPI_Type_commit(&bcInfoType);
	MPI_Bcast(bc[0], 6, bcInfoType, 0, PETSC_COMM_WORLD);
	MPI_Bcast(bc[1], 6, bcInfoType, 0, PETSC_COMM_WORLD);
	MPI_Bcast(bc[2], 6, bcInfoType, 0, PETSC_COMM_WORLD);
	MPI_Type_free(&bcInfoType);
}
