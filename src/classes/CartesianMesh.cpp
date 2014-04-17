#include <fstream>
#include "CartesianMesh.h"
#include "yaml-cpp/yaml.h"

template <PetscInt dim>
CartesianMesh<dim>::CartesianMesh(std::string fileName)
{
	PetscInt       rank;
	
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	// first pass
	if(rank==0)
	{
		PetscReal     start;
		PetscInt      numCells;
		std::string   direction;
		std::ifstream file(fileName.c_str());
		YAML::Parser  parser(file);
		YAML::Node    doc;
		
		parser.GetNextDocument(doc);
				
		// cycle through each direction
		for (size_t i=0; i<doc.size(); i++)
		{			
			doc[i]["direction"] >> direction;
			doc[i]["start"] >> start;
				
			if(direction == "x")
				nx = 0;
			else if(direction == "y")
				ny = 0;
			else if(direction == "z")
				nz = 0;
			
			const YAML::Node &subDomains = doc[i]["subDomains"];
			for (size_t j=0; j<subDomains.size(); j++)
			{
				subDomains[j]["cells"] >> numCells;
				if (direction == "x")
					nx += numCells;
				else if(direction == "y")
					ny += numCells;
				else if(direction == "z")
					nz += numCells;
			}
		}
	}
	MPI_Barrier(PETSC_COMM_WORLD);
	
	// broadcast vector sizes to all processes
	MPI_Bcast(&nx, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
	MPI_Bcast(&ny, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
	MPI_Bcast(&nz, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
	
	// allocate memory
	x.resize(nx+1);
	dx.resize(nx);
	
	y.resize(ny+1);
	dy.resize(ny);
	
	z.resize(nz+1);	
	dz.resize(nz);
	
	// second pass
	if(rank == 0)
	{
		size_t        numCells, first;
		std::string   direction;
		std::ifstream file(fileName.c_str());
		YAML::Parser  parser(file);
		YAML::Node    doc;
		PetscReal     start,
		              end,
			          stretchRatio,
			          h;
		
		parser.GetNextDocument(doc);
				
		// cycle through each direction
		for (size_t i=0; i<doc.size(); i++)
		{
			doc[i]["direction"] >> direction;
			doc[i]["start"] >> start;
			
			const YAML::Node &subDomains = doc[i]["subDomains"];
			
			first = 0;
			if(direction=="x")
				x[first] = start;
			else if(direction=="y")
				y[first] = start;
			else if(direction=="z")
				z[first] = start;
			
			for (size_t i=0; i<subDomains.size(); i++)
			{
				subDomains[i]["end"] >> end;
				subDomains[i]["cells"] >> numCells;
				subDomains[i]["stretchRatio"] >> stretchRatio;
				
				if(fabs(stretchRatio-1.0) < 1.0e-6)
				{
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
					}
				}
				else
				{
					h = (end - start)*(stretchRatio-1)/(pow(stretchRatio, numCells)-1);
					for(size_t j=first; j<first+numCells; j++)
					{
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
					}
				}
				first += numCells;
				start = end;
			}
		}
	}
	MPI_Barrier(PETSC_COMM_WORLD);
	
	// broadcast vectors to all processes
	MPI_Bcast(&x.front(), nx+1, MPIU_REAL, 0, MPI_COMM_WORLD);
	MPI_Bcast(&dx.front(), nx, MPIU_REAL, 0, MPI_COMM_WORLD);
	
	MPI_Bcast(&y.front(), ny+1, MPIU_REAL, 0, MPI_COMM_WORLD);
	MPI_Bcast(&dy.front(), ny, MPIU_REAL, 0, MPI_COMM_WORLD);
	
	MPI_Bcast(&z.front(), nz+1, MPIU_REAL, 0, MPI_COMM_WORLD);	
	MPI_Bcast(&dz.front(), nz, MPIU_REAL, 0, MPI_COMM_WORLD);
	
	MPI_Barrier(PETSC_COMM_WORLD);
	
	std::cout << "nx: " << nx << ", ny: " << ny << std::endl;
	for(std::vector<PetscReal>::iterator it=x.begin(); it!=x.end(); ++it)
	{
		std::cout << *it << ", ";
	}
	std::cout << std::endl;
}

template class CartesianMesh<2>;
template class CartesianMesh<3>;
