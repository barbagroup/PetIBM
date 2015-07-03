template <>
PetscErrorCode TairaColoniusSolver<2>::writeForces()
{
	PetscErrorCode ierr;
	PetscInt       rank;

	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

	if(rank==0)
	{
		std::string filename = caseFolder + "/forces.txt";
		if(timeStep==1)
		{
			forcesFile.open(filename.c_str());
		}
		else
		{
			forcesFile.open(filename.c_str(), std::ios::out | std::ios::app);
		}
		forcesFile << timeStep*simParams->dt << '\t' << force[0] << '\t' << force[1] << std::endl;
		forcesFile.close();
	}

	return 0;
} // writeForces

template <>
PetscErrorCode TairaColoniusSolver<3>::writeForces()
{
	PetscErrorCode ierr;
	PetscInt       rank;

	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

	if(rank==0)
	{
		std::string filename = caseFolder + "/forces.txt";
		if(timeStep==1)
		{
			forcesFile.open(filename.c_str());
		}
		else
		{
			forcesFile.open(filename.c_str(), std::ios::out | std::ios::app);
		}
		forcesFile << timeStep*simParams->dt << '\t' << force[0] << '\t' << force[1] << '\t' << force[2] << std::endl;
		forcesFile.close();
	}

	return 0;
} // writeForces
