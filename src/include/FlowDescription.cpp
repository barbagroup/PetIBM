/***************************************************************************//**
* \file
* \brief Source file to define member functions of FlowDescription
*/

#include "FlowDescription.h"
#include "yaml-cpp/yaml.h"
#include <iostream>
#include <fstream>

/***************************************************************************//**
* \brief Converts \c std::string to \c Boundary
*/
Boundary boundaryFromString(std::string &s)
{
	if (s == "xMinus") return XMINUS;
	if (s == "xPlus")  return XPLUS;
	if (s == "yMinus") return YMINUS;
	if (s == "yPlus")  return YPLUS;
	if (s == "zMinus") return ZMINUS;
	if (s == "zPlus")  return ZPLUS;
	
	std::cout << "ERROR: Invalid boundary location!\n";
	exit(0);
}

/***************************************************************************//**
* \brief Converts \c std::string to \c BCType
*/
BCType bcTypeFromString(std::string &s)
{
	if (s == "DIRICHLET") return DIRICHLET;
	if (s == "NEUMANN") return NEUMANN;
	if (s == "CONVECTIVE") return CONVECTIVE;
	if (s == "PERIODIC") return PERIODIC;
	
	std::cout << "ERROR: Invalid boundary condition type!\n";
	exit(0);
}

/***************************************************************************//**
* \param fileName Input file path
*
* This is the constructor for the class FlowDescription. A case folder with
* the input files is supplied to the flow solver, and this function reads the
* file \c flowDescription.yaml in the folder. The parameters are listed in
* the file using the YAML format.
*/
FlowDescription::FlowDescription(std::string fileName)
{
	PetscInt    rank;
	
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank); // get rank of current process
	
	if(rank==0) // read the input file only on process 0
	{
		std::ifstream inputFile(fileName.c_str());
		YAML::Parser  parser(inputFile);
		YAML::Node    document;
		std::string   locationString, // store the boundary location (+/-X, +/-Y, +/-Z) in a string
		              typeString;     // stores the type of boundary condition in a string
		Boundary      location;
		
		parser.GetNextDocument(document);
		
		document[0]["dimensions"] >> dimensions; // read the number of dimensions
		document[0]["nu"] >> nu;                 // read the kinematic viscosity
		
		document[0]["initialVelocity"][0] >> initialVelocity[0]; // read the x-component of initial velocity
		document[0]["initialVelocity"][1] >> initialVelocity[1]; // read the y-component of initial velocity
		if(dimensions==3) // flow 3-D flows
		{
			document[0]["initialVelocity"][2] >> initialVelocity[2]; // read the z-component of initial velocity
		}
		else
		{
			initialVelocity[2] = 0.0; // set the z-component of initial velocity to 0 for 2-D flows
		}
		
		// specifying the initial perturbation is optional
		// hence it is inside a try block
		try
		{
			document[0]["initialPerturbation"][0] >> initialPerturbation[0]; // read the initial perturbation in u
			document[0]["initialPerturbation"][1] >> initialPerturbation[1]; // read the initial perturbation in v
			if(dimensions==3) // for 3-D flows
			{
				document[0]["initialPerturbation"][2] >> initialPerturbation[2]; // read the initial perturbation in w
			}
			else
			{
				initialVelocity[2] = 0.0; // set the initial perturbation in w to 0 for 2-D flows
			}
		}
		catch(...)
		{
		}
		
		const YAML::Node &bcNode = document[0]["boundaryConditions"];
		
		for (size_t i=0; i<bcNode.size(); i++) // cycle through all boundaries
		{
			bcNode[i]["location"] >> locationString;             // read which boundary is being described
			location = boundaryFromString(locationString);       // set the boundary by converting from string
			
			bcNode[i]["u"][0] >> typeString;                     // read the type of boundary condition for u
			bc[0][location].type = bcTypeFromString(typeString); // set the type by converting from string
			bcNode[i]["u"][1] >> bc[0][location].value;          // read the value associated with it
			
			bcNode[i]["v"][0] >> typeString;                     // read the type of boundary condition for v
			bc[1][location].type = bcTypeFromString(typeString); // set the type by converting from string
			bcNode[i]["v"][1] >> bc[1][location].value;          // read the value associated with it
			
			if(dimensions==3) // for 3-D flows
			{
				bcNode[i]["w"][0] >> typeString;                     // read the type of boundary condition for v
				bc[2][location].type = bcTypeFromString(typeString); // set the type by converting from string
				bcNode[i]["w"][1] >> bc[2][location].value;          // read the value associated with it
			}
			else
			{
				bc[2][location].type  = PERIODIC; // set the type to periodic for 2-D flows
				bc[2][location].value = 0.0;
			}
		}
		
		// run some sanity checks on the input data
		PetscBool flag = PETSC_TRUE;
		// if the boundary condition on one face is periodic,
		// then it should be periodic on the opposite face too
		// check u on the X and Y faces
		if(bc[0][XMINUS].type==PERIODIC && bc[0][XPLUS].type!=PERIODIC) flag = PETSC_FALSE;
		if(bc[0][XMINUS].type!=PERIODIC && bc[0][XPLUS].type==PERIODIC) flag = PETSC_FALSE;
		if(bc[0][YMINUS].type==PERIODIC && bc[0][YPLUS].type!=PERIODIC) flag = PETSC_FALSE;
		if(bc[0][YMINUS].type!=PERIODIC && bc[0][YPLUS].type==PERIODIC) flag = PETSC_FALSE;
		if(dimensions==3)
		{
			// check u on the Z faces
			if(bc[0][ZMINUS].type==PERIODIC && bc[0][ZPLUS].type!=PERIODIC) flag = PETSC_FALSE;
			if(bc[0][ZMINUS].type!=PERIODIC && bc[0][ZPLUS].type==PERIODIC) flag = PETSC_FALSE;
		}
		// if the boundary condition for one component of velocity is periodic,
		// then it should be periodic for the other components too
		// compare u and v on the X and Y faces
		if(bc[0][XMINUS].type==PERIODIC && bc[1][XMINUS].type!=PERIODIC) flag = PETSC_FALSE;
		if(bc[0][XPLUS].type==PERIODIC && bc[1][XPLUS].type!=PERIODIC) flag = PETSC_FALSE;
		if(bc[0][YMINUS].type==PERIODIC && bc[1][YMINUS].type!=PERIODIC) flag = PETSC_FALSE;
		if(bc[0][YPLUS].type==PERIODIC && bc[1][YPLUS].type!=PERIODIC) flag = PETSC_FALSE;
		if(dimensions==3)
		{
			// compare u and w on the X and Y faces
			if(bc[0][XMINUS].type==PERIODIC && bc[2][XMINUS].type!=PERIODIC) flag = PETSC_FALSE;
			if(bc[0][XPLUS].type==PERIODIC && bc[2][XPLUS].type!=PERIODIC) flag = PETSC_FALSE;
			if(bc[0][YMINUS].type==PERIODIC && bc[2][YMINUS].type!=PERIODIC) flag = PETSC_FALSE;
			if(bc[0][YPLUS].type==PERIODIC && bc[2][YPLUS].type!=PERIODIC) flag = PETSC_FALSE;
			
			// compare u and v, u and w the on Z faces
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
		inputFile.close();
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
