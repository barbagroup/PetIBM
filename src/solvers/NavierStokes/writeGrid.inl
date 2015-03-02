template <>
PetscErrorCode NavierStokesSolver<2>::writeGrid()
{
	PetscErrorCode ierr;
	PetscInt       rank;
	
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
	
	if(rank==0)
	{
		std::ofstream f(caseFolder+"/grid.txt");
		f << mesh->nx << '\t' << mesh->ny << '\n';
		for(std::vector<PetscReal>::const_iterator i=mesh->x.begin(); i!=mesh->x.end(); ++i)
			f << *i << '\n';
		for(std::vector<PetscReal>::const_iterator i=mesh->y.begin(); i!=mesh->y.end(); ++i)
			f << *i << '\n';
		f.close();
	}

	return 0;
}

template <>
PetscErrorCode NavierStokesSolver<3>::writeGrid()
{
	PetscErrorCode ierr;
	PetscInt       rank;
	
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
	
	if(rank==0)
	{
		std::ofstream f(caseFolder+"/grid.txt");
		f << mesh->nx << '\t' << mesh->ny << '\t' << mesh->nz << '\n';
		for(std::vector<PetscReal>::const_iterator i=mesh->x.begin(); i!=mesh->x.end(); ++i)
			f << *i << '\n';
		for(std::vector<PetscReal>::const_iterator i=mesh->y.begin(); i!=mesh->y.end(); ++i)
			f << *i << '\n';
		for(std::vector<PetscReal>::const_iterator i=mesh->z.begin(); i!=mesh->z.end(); ++i)
			f << *i << '\n';
		f.close();
	}

	return 0;
}