template <>
void NavierStokesSolver<2>::generateR2()
{
	PetscErrorCode ierr;
	PetscInt       mstart, nstart, m, n, i, j;
	PetscReal      **qx, **qy;
	PetscReal      **bc2;
	Vec            bc2Global;

	ierr = VecSet(r2, 0.0); CHKERRV(ierr);
	ierr = DMCompositeGetAccess(phiPack, r2,  &bc2Global); CHKERRV(ierr);

	ierr = DMDAVecGetArray(uda, qxLocal, &qx); CHKERRV(ierr);
	ierr = DMDAVecGetArray(vda, qyLocal, &qy); CHKERRV(ierr);

	ierr = DMDAVecGetArray(pda, bc2Global, &bc2); CHKERRV(ierr);
	ierr = DMDAGetCorners(pda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRV(ierr);
	// x-faces
	if(flowDesc->bc[0][XPLUS].type != PERIODIC) // don't update if the BC type is periodic
	{
		//loop over all points on the x-face
		for(j=nstart; j<nstart+n; j++)
		{
			// -X
			if(mstart == 0) // if the -X face is in the current process
			{
				bc2[j][0] -= qx[j][-1];
			}
			// +X
			if(mstart+m-1 == mesh->nx-1) // if the +X face is in the current process
			{
				bc2[j][mesh->nx-1] += qx[j][mesh->nx-1];
			}
		}
	}
	// y-faces
	if(flowDesc->bc[1][YPLUS].type != PERIODIC) // don't update if the BC type is periodic
	{
		//loop over all points on the y-face
		for(i=mstart; i<mstart+m; i++)
		{
			// -Y
			if(nstart == 0) // if the -Y face is in the current process
			{
				bc2[0][i] -= qy[-1][i];
			}
			// +Y
			if(nstart+n-1 == mesh->ny-1) // if the +Y face is in the current process
			{
				bc2[mesh->ny-1][i] += qy[mesh->ny-1][i];
			}
		}
	}
	ierr = DMDAVecRestoreArray(pda, bc2Global, &bc2); CHKERRV(ierr);

	ierr = DMDAVecRestoreArray(uda, qxLocal, &qx); CHKERRV(ierr);
	ierr = DMDAVecRestoreArray(vda, qyLocal, &qy); CHKERRV(ierr);

	ierr = DMCompositeRestoreAccess(phiPack, r2,  &bc2Global); CHKERRV(ierr);
}

template <>
void NavierStokesSolver<3>::generateR2()
{
}