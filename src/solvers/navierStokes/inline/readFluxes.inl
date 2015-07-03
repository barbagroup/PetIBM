template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::readFluxes()
{
	PetscErrorCode    ierr;
	PetscViewer       viewer;
	std::stringstream ss;
	std::string       savePointDir, fileName;
	Vec               qxGlobal, qyGlobal, qzGlobal;

	// get access to the individual vectors of the composite vector
	// depending on whether it is a 2-D or a 3-D flow
	if(dim==2)
	{
		ierr = DMCompositeGetAccess(qPack, q, &qxGlobal, &qyGlobal); CHKERRQ(ierr);
	}
	else if(dim==3)
	{
		ierr = DMCompositeGetAccess(qPack, q, &qxGlobal, &qyGlobal, &qzGlobal); CHKERRQ(ierr);
	}

	ierr = PetscPrintf(PETSC_COMM_WORLD, "\nRestarting from time step %d.\n", timeStep); CHKERRQ(ierr);

	// the name of the folder is the time step at which data is saved
	// 7 characters long, with leading zeros
	ss << caseFolder << "/" << std::setfill('0') << std::setw(7) << timeStep;
	savePointDir = ss.str();

	// save x-component of fluxes
	fileName = savePointDir + "/qx.dat";
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, fileName.c_str(), FILE_MODE_READ, &viewer); CHKERRQ(ierr);
	ierr = VecLoad(qxGlobal, viewer); CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

	// save y-component of fluxes
	fileName = savePointDir + "/qy.dat";
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, fileName.c_str(), FILE_MODE_READ, &viewer); CHKERRQ(ierr);
	ierr = VecLoad(qyGlobal, viewer); CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

	// save z-component of fluxes if it is a 3-D flow
	if(dim==3)
	{
		fileName = savePointDir + "/qz.dat";
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, fileName.c_str(), FILE_MODE_READ, &viewer); CHKERRQ(ierr);
		ierr = VecLoad(qzGlobal, viewer); CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
	}

	if(dim==2)
	{
		ierr = DMCompositeRestoreAccess(qPack, q, &qxGlobal, &qyGlobal); CHKERRQ(ierr);
	}
	else if(dim==3)
	{
		ierr = DMCompositeRestoreAccess(qPack, q, &qxGlobal, &qyGlobal, &qzGlobal); CHKERRQ(ierr);
	}
	
	return 0;
} // readFluxes
