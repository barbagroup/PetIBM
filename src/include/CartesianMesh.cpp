/***************************************************************************//**
* \file
* \brief Source file to define member functions of CartesianMesh
*/

#include "CartesianMesh.h"
#include "yaml-cpp/yaml.h"
#include <fstream>

CartesianMesh::CartesianMesh()
{
}

CartesianMesh::CartesianMesh(std::string fileName)
{
	initialize(fileName);
}

/***************************************************************************//**
* \param fileName Input file path
*
* This function initializes an object of type CartesianMesh. A case folder with
* the input files is supplied to the flow solver, and this function reads the
* file \c cartesianMesh.yaml in the folder. The mesh is described in the file
* using the YAML format.
*/
void CartesianMesh::initialize(std::string fileName)
{
	PetscInt rank;
	
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank); // get rank of the current process
	
	nx = 0;
	ny = 0;
	nz = 0;
	
	// first pass of the input file
	// do this to figure out the number of cells in each direction
	if(rank==0)
	{
		PetscInt      numCells;
		std::string   direction;
		std::ifstream inputFile(fileName.c_str());
		YAML::Parser  parser(inputFile);
		YAML::Node    document;
		
		parser.GetNextDocument(document);
				
		// cycle through each direction
		for (size_t i=0; i<document.size(); i++)
		{			
			document[i]["direction"] >> direction; // read the direction ("x", "y" or "z")

			// initialize the number of cells to zero
			if(direction == "x") nx = 0;
			else if(direction == "y") ny = 0;
			else if(direction == "z") nz = 0;
			
			// read the number of cells in each subdomain
			// and add it to the total number of cells
			const YAML::Node &subDomains = document[i]["subDomains"];
			for (size_t j=0; j<subDomains.size(); j++)
			{
				subDomains[j]["cells"] >> numCells;
				if (direction == "x") nx += numCells;
				else if(direction == "y") ny += numCells;
				else if(direction == "z") nz += numCells;
			}
		}
	}
	MPI_Barrier(PETSC_COMM_WORLD);
	
	// broadcast number of cells to all processes
	MPI_Bcast(&nx, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
	MPI_Bcast(&ny, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
	MPI_Bcast(&nz, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
	
	// allocate memory
	// number of nodes will be 1 greater than the number of cells
	x.resize(nx+1);
	dx.resize(nx);
	y.resize(ny+1);
	dy.resize(ny);
	if(nz > 0)
	{
		z.resize(nz+1);
		dz.resize(nz);
	}
	
	// second pass of the input file
	// to calculate the coordinates of the nodes, and the cell widths
	if(rank == 0)
	{
		size_t        numCells,
		              first;
		std::string   direction;
		std::ifstream inputFile(fileName.c_str());
		YAML::Parser  parser(inputFile);
		YAML::Node    document;
		PetscReal     start,
		              end,
			          stretchRatio,
			          h;
		
		parser.GetNextDocument(document);
				
		// cycle through each direction
		for (size_t i=0; i<document.size(); i++)
		{
			document[i]["direction"] >> direction; // read the direction
			document[i]["start"] >> start;         // coordinate of the first node in the subdomain
			
			const YAML::Node &subDomains = document[i]["subDomains"];
			
			// set the coordinate of the first node in the first subdomain
			first = 0;
			if(direction=="x") x[first] = start;
			else if(direction=="y") y[first] = start;
			else if(direction=="z") z[first] = start;
			
			// cycle through the subdomains
			for (size_t i=0; i<subDomains.size(); i++)
			{
				subDomains[i]["end"] >> end; // read the coordinate of the last node in the subdomain
				subDomains[i]["cells"] >> numCells; // read the number of cells in the subdomain
				subDomains[i]["stretchRatio"] >> stretchRatio; // read the stretching ratio of the cells in subdomain
				
				if(fabs(stretchRatio-1.0) < 1.0e-6) // if the stretching ratio is 1
				{
					// the cells are of uniform width h
					h = (end - start)/numCells;
					for(size_t j=first; j<first+numCells; j++)
					{
						if(direction=="x")
						{
							dx[j]  = h;
							x[j+1] = x[j] + dx[j];
						}
						else if(direction=="y")
						{
							dy[j]  = h;
							y[j+1] = y[j] + dy[j];
						}
						else if(direction=="z")
						{
							dz[j]  = h;
							z[j+1] = z[j] + dz[j];
						}
					}
				}
				else // if the stretching ratio is different from 1
				{
					// width of the first cell in the subdomain
					h = (end - start)*(stretchRatio-1)/(pow(stretchRatio, numCells)-1);
					for(size_t j=first; j<first+numCells; j++)
					{
						// obtain the widths of subsequent cells by
						// multiplying by the stretching ratio
						if(direction=="x")
						{
							dx[j]  = h*pow(stretchRatio, j-first);
							x[j+1] = x[j] + dx[j];
						}
						else if(direction=="y")
						{
							dy[j]  = h*pow(stretchRatio, j-first);
							y[j+1] = y[j] + dy[j];
						}
						else if(direction=="z")
						{
							dz[j]  = h*pow(stretchRatio, j-first);
							z[j+1] = z[j] + dz[j];
						}
					}
				}
				// the first node in the next subdomain is the last node in the current subdomain
				first += numCells;
				start = end;
			}
		}
	}
	MPI_Barrier(PETSC_COMM_WORLD);
	
	// broadcast vectors to all processes
	MPI_Bcast(&x.front(), nx+1, MPIU_REAL, 0, PETSC_COMM_WORLD);
	MPI_Bcast(&dx.front(), nx, MPIU_REAL, 0, PETSC_COMM_WORLD);
	MPI_Bcast(&y.front(), ny+1, MPIU_REAL, 0, PETSC_COMM_WORLD);
	MPI_Bcast(&dy.front(), ny, MPIU_REAL, 0, PETSC_COMM_WORLD);
	if(nz > 0)
	{
		MPI_Bcast(&z.front(), nz+1, MPIU_REAL, 0, PETSC_COMM_WORLD);
		MPI_Bcast(&dz.front(), nz, MPIU_REAL, 0, PETSC_COMM_WORLD);
	}
}
