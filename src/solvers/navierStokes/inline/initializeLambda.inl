template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::initializeLambda()
{
	PetscErrorCode ierr;
	Vec            phi;
	
	ierr = DMCompositeGetAccess(lambdaPack, lambda, &phi); CHKERRQ(ierr);

	if (simParams->startStep > 0 || simParams->restartFromSolution)
	{
		PetscViewer       viewer;
		std::stringstream ss;
		std::string       savePointDir, fileName;

		// the name of the folder is the time step at which data is saved
		// 7 characters long, with leading zeros
		ss << caseFolder << "/" << std::setfill('0') << std::setw(7) << timeStep;
		savePointDir = ss.str();

		fileName = savePointDir + "/phi.dat";
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, fileName.c_str(), FILE_MODE_READ, &viewer); CHKERRQ(ierr);
		ierr = VecLoad(phi, viewer); CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
	}
	
	ierr = DMCompositeRestoreAccess(lambdaPack, lambda, &phi); CHKERRQ(ierr);

	return 0;
} // initializeLambda