template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::initializeLambda()
{
	PetscErrorCode ierr;
	Vec            phi, fTilde;
	
	ierr = DMCompositeGetAccess(NavierStokesSolver<dim>::lambdaPack, NavierStokesSolver<dim>::lambda, &phi, &fTilde); CHKERRQ(ierr);

	if(NavierStokesSolver<dim>::simParams->restart)
	{
		PetscViewer       viewer;
		std::stringstream ss;
		std::string       savePointDir, fileName;

		ss << NavierStokesSolver<dim>::caseFolder << "/" << std::setfill('0') << std::setw(7) << NavierStokesSolver<dim>::timeStep;
		savePointDir = ss.str();

		ss.str("");
		ss.clear();
		ss << savePointDir << "/phi.dat";
		fileName = ss.str();
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, fileName.c_str(), FILE_MODE_READ, &viewer); CHKERRQ(ierr);
		ierr = VecLoad(phi, viewer); CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

		ss.str("");
		ss.clear();
		ss << savePointDir << "/fTilde.dat";
		fileName = ss.str();
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, fileName.c_str(), FILE_MODE_READ, &viewer); CHKERRQ(ierr);
		ierr = VecLoad(fTilde, viewer); CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
	}
	
	ierr = DMCompositeRestoreAccess(NavierStokesSolver<dim>::lambdaPack, NavierStokesSolver<dim>::lambda, &phi, &fTilde); CHKERRQ(ierr);

	return 0;
}