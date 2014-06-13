#include <sstream>

template <>
PetscErrorCode NavierStokesSolver<2>::initialiseFluxes()
{
	PetscErrorCode ierr;
	Vec            qxGlobal, qyGlobal;
	
	ierr = DMCompositeGetAccess(qPack, q, &qxGlobal, &qyGlobal); CHKERRQ(ierr);

	if(simParams->restart)
	{
		PetscViewer       viewer;
		std::stringstream ss;
		std::string       savePointDir, fileName;

		PetscPrintf(PETSC_COMM_WORLD, "Restarting from time step %d.\n", timeStep);

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
	}
	else
	{
		PetscInt       mstart, nstart, m, n;
		PetscReal      **qx, **qy;

		// U-FLUXES
		ierr = DMDAVecGetArray(uda, qxGlobal, &qx); CHKERRQ(ierr);
		ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
		// Set interior values for u-fluxes
		for(PetscInt j=nstart; j<nstart+n; j++)
		{
			for(PetscInt i=mstart; i<mstart+m; i++)
			{
				qx[j][i] = flowDesc->initialVelocity[0] * mesh->dy[j];
			}
		}
		ierr = DMDAVecRestoreArray(uda, qxGlobal, &qx); CHKERRQ(ierr);
		
		// V-FLUXES
		ierr = DMDAVecGetArray(vda, qyGlobal, &qy); CHKERRQ(ierr);
		ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
		// Set interior values for v-fluxes
		for(PetscInt j=nstart; j<nstart+n; j++)
		{
			for(PetscInt i=mstart; i<mstart+m; i++)
			{
				qy[j][i] = flowDesc->initialVelocity[1] * mesh->dx[i];
			}
		}
		ierr = DMDAVecRestoreArray(vda, qyGlobal, &qy); CHKERRQ(ierr);
	}
	
	ierr = DMCompositeRestoreAccess(qPack, q, &qxGlobal, &qyGlobal); CHKERRQ(ierr);
	ierr = DMCompositeScatter(qPack, q, qxLocal, qyLocal); CHKERRQ(ierr);

	return 0;
}

template <>
PetscErrorCode NavierStokesSolver<3>::initialiseFluxes()
{
	PetscErrorCode ierr;
	Vec            qxGlobal, qyGlobal, qzGlobal;
	
	ierr = DMCompositeGetAccess(qPack, q, &qxGlobal, &qyGlobal, &qzGlobal); CHKERRQ(ierr);

	if(simParams->restart)
	{
		PetscViewer       viewer;
		std::stringstream ss;
		std::string       savePointDir, fileName;

		PetscPrintf(PETSC_COMM_WORLD, "Restarting from time step %d.\n", timeStep);

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

		ss.str("");
		ss.clear();
		ss << savePointDir << "/qz.dat";
		fileName = ss.str();
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, fileName.c_str(), FILE_MODE_READ, &viewer); CHKERRQ(ierr);
		ierr = VecLoad(qzGlobal, viewer); CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
	}
	else
	{
		PetscInt       mstart, nstart, pstart, m, n, p;
		PetscReal      ***qx, ***qy, ***qz;
		
		// U-FLUXES
		ierr = DMDAVecGetArray(uda, qxGlobal, &qx); CHKERRQ(ierr);
		ierr = DMDAGetCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
		// Set interior values for u-fluxes
		for(PetscInt k=pstart; k<pstart+p; k++)
		{
			for(PetscInt j=nstart; j<nstart+n; j++)
			{
				for(PetscInt i=mstart; i<mstart+m; i++)
				{
					qx[k][j][i] = flowDesc->initialVelocity[0] * (mesh->dy[j]*mesh->dz[k]);
				}
			}
		}
		ierr = DMDAVecRestoreArray(uda, qxGlobal, &qx); CHKERRQ(ierr);
		
		// V-FLUXES
		ierr = DMDAVecGetArray(vda, qyGlobal, &qy); CHKERRQ(ierr);
		ierr = DMDAGetCorners(vda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
		// Set interior values for v-fluxes
		for(PetscInt k=pstart; k<pstart+p; k++)
		{
			for(PetscInt j=nstart; j<nstart+n; j++)
			{
				for(PetscInt i=mstart; i<mstart+m; i++)
				{
					qy[k][j][i] = flowDesc->initialVelocity[1] * (mesh->dx[i]*mesh->dz[k]);
				}
			}
		}
		ierr = DMDAVecRestoreArray(vda, qyGlobal, &qy); CHKERRQ(ierr);

		// W-FLUXES
		ierr = DMDAVecGetArray(wda, qzGlobal, &qz); CHKERRQ(ierr);
		ierr = DMDAGetCorners(wda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
		// Set interior values for w-fluxes
		for(PetscInt k=pstart; k<pstart+p; k++)
		{
			for(PetscInt j=nstart; j<nstart+n; j++)
			{
				for(PetscInt i=mstart; i<mstart+m; i++)
				{
					qz[k][j][i] = flowDesc->initialVelocity[2] * (mesh->dx[i]*mesh->dy[j]);
				}
			}
		}
		ierr = DMDAVecRestoreArray(wda, qzGlobal, &qz); CHKERRQ(ierr);
	}
	
	ierr = DMCompositeRestoreAccess(qPack, q, &qxGlobal, &qyGlobal, &qzGlobal); CHKERRQ(ierr);
	ierr = DMCompositeScatter(qPack, q, qxLocal, qyLocal, qzLocal); CHKERRQ(ierr);

	return 0;
}
