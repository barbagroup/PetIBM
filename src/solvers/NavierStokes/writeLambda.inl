template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::writeLambda()
{
	PetscErrorCode  ierr;
	Vec             phi;
	std::string     savePointDir, fileName;
	PetscViewer     viewer;
	
	// create output folder
	std::stringstream ss;
	ss << caseFolder << "/" << std::setfill('0') << std::setw(7) << timeStep;
	savePointDir = ss.str();

	ierr = DMCompositeGetAccess(lambdaPack, lambda, &phi); CHKERRQ(ierr);

	// print phi to file
	ss.str("");
	ss.clear();
	ss << savePointDir << "/phi.dat";
	fileName = ss.str();
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, fileName.c_str(), FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
	ierr = VecView(phi, viewer); CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

	ierr = DMCompositeRestoreAccess(lambdaPack, lambda, &phi); CHKERRQ(ierr);

	return 0;
}