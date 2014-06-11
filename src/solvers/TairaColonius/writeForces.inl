template <>
void TairaColoniusSolver<2>::writeForces()
{
	PetscInt        rank;

	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	if(rank==0)
	{
		std::string filename = caseFolder + "/forces.txt";
		if(timeStep==1)
			forcesFile.open(filename.c_str());
		else	
			forcesFile.open(filename.c_str(), std::ios::out | std::ios::app);
		forcesFile << timeStep*simParams->dt << '\t' << force[0] << '\t' << force[1] << std::endl;
		forcesFile.close();
	}
}

template <>
void TairaColoniusSolver<3>::writeForces()
{
	PetscInt        rank;

	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	if(rank==0)
	{
		std::string filename = caseFolder + "/forces.txt";
		if(timeStep==1)
			forcesFile.open(filename.c_str());
		else	
			forcesFile.open(filename.c_str(), std::ios::out | std::ios::app);
		forcesFile << timeStep*simParams->dt << '\t' << force[0] << '\t' << force[1] << '\t' << force[2] << std::endl;
		forcesFile.close();
	}
}