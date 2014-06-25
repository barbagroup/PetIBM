template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::initializeLambda()
{
	PetscErrorCode ierr;
	Vec            phi;
	
	ierr = DMCompositeGetAccess(lambdaPack, lambda, &phi); CHKERRQ(ierr);

	if(simParams->restart)
	{
		PetscViewer       viewer;
		std::stringstream ss;
		std::string       savePointDir, fileName;

		ss << caseFolder << "/" << std::setfill('0') << std::setw(7) << timeStep;
		savePointDir = ss.str();

		ss.str("");
		ss.clear();
		ss << savePointDir << "/phi.dat";
		fileName = ss.str();
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, fileName.c_str(), FILE_MODE_READ, &viewer); CHKERRQ(ierr);
		ierr = VecLoad(phi, viewer); CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
	}
	
	ierr = DMCompositeRestoreAccess(lambdaPack, lambda, &phi); CHKERRQ(ierr);

	return 0;
}