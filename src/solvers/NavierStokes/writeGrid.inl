#include <fstream>

template <>
void NavierStokesSolver<2>::writeGrid()
{
	PetscInt    rank;
	
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	if(rank==0)
	{
		std::ofstream f(caseFolder+"/grid.txt");
		for(std::vector<PetscReal>::const_iterator i=mesh->x.begin(); i!=mesh->x.end(); ++i)
			f << *i << '\n';
		for(std::vector<PetscReal>::const_iterator i=mesh->y.begin(); i!=mesh->y.end(); ++i)
			f << *i << '\n';
		f.close();
	}
}

template <>
void NavierStokesSolver<3>::writeGrid()
{
	PetscInt    rank;
	
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	if(rank==0)
	{
		std::ofstream f(caseFolder+"/grid.txt");
		for(std::vector<PetscReal>::const_iterator i=mesh->x.begin(); i!=mesh->x.end(); ++i)
			f << *i << '\n';
		for(std::vector<PetscReal>::const_iterator i=mesh->y.begin(); i!=mesh->y.end(); ++i)
			f << *i << '\n';
		for(std::vector<PetscReal>::const_iterator i=mesh->z.begin(); i!=mesh->z.end(); ++i)
			f << *i << '\n';
		f.close();
	}
}