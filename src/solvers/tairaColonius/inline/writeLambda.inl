template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::writeLambda()
{
	PetscErrorCode  ierr;
	Vec             phi, fTilde;
	std::string     savePointDir, fileName;
	PetscViewer     viewer;
	
	// create output folder
	std::stringstream ss;
	ss << NavierStokesSolver<dim>::caseFolder << "/" << std::setfill('0') << std::setw(7) << NavierStokesSolver<dim>::timeStep;
	savePointDir = ss.str();

	ierr = DMCompositeGetAccess(NavierStokesSolver<dim>::lambdaPack, NavierStokesSolver<dim>::lambda, &phi, &fTilde); CHKERRQ(ierr);

	// print phi to file
	ss.str("");
	ss.clear();
	ss << savePointDir << "/phi.dat";
	fileName = ss.str();
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, fileName.c_str(), FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
	ierr = VecView(phi, viewer); CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

	// print fTilde to file
	ss.str("");
	ss.clear();
	ss << savePointDir << "/fTilde.dat";
	fileName = ss.str();
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, fileName.c_str(), FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
	ierr = VecView(fTilde, viewer); CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

	ierr = DMCompositeRestoreAccess(NavierStokesSolver<dim>::lambdaPack, NavierStokesSolver<dim>::lambda, &phi, &fTilde); CHKERRQ(ierr);

	return 0;
} // writeLambda