template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::readFluxes(Vec qxGlobal, Vec qyGlobal, Vec qzGlobal)
{
	PetscErrorCode    ierr;
	PetscViewer       viewer;
	std::stringstream ss;
	std::string       savePointDir, fileName;

	ierr = PetscPrintf(PETSC_COMM_WORLD, "Restarting from time step %d.\n", timeStep); CHKERRQ(ierr);

	ss << caseFolder << "/" << std::setfill('0') << std::setw(7) << timeStep;
	savePointDir = ss.str();

	ss.str("");
	ss.clear();
	ss << savePointDir << "/qx.dat";
	fileName = ss.str();
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, fileName.c_str(), FILE_MODE_READ, &viewer); CHKERRQ(ierr);
	ierr = VecLoad(qxGlobal, viewer); CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

	ss.str("");
	ss.clear();
	ss << savePointDir << "/qy.dat";
	fileName = ss.str();
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, fileName.c_str(), FILE_MODE_READ, &viewer); CHKERRQ(ierr);
	ierr = VecLoad(qyGlobal, viewer); CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

	if(dim==3)
	{
		ss.str("");
		ss.clear();
		ss << savePointDir << "/qz.dat";
		fileName = ss.str();
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, fileName.c_str(), FILE_MODE_READ, &viewer); CHKERRQ(ierr);
		ierr = VecLoad(qzGlobal, viewer); CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
	}
	
	return 0;
}