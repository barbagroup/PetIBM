template <>
void NavierStokesSolver<2>::generateR2()
{
	PetscErrorCode ierr;
	PetscInt       mstart, nstart, m, n, i, j;
	PetscReal      **qx, **qy;
	PetscReal      **bc2;
	Vec            bc2Global;

	ierr = VecSet(r2, 0.0); CHKERRV(ierr);
	ierr = DMCompositeGetAccess(lambdaPack, r2,  &bc2Global); CHKERRV(ierr);

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

	ierr = DMCompositeRestoreAccess(lambdaPack, r2,  &bc2Global); CHKERRV(ierr);
}

template <>
void NavierStokesSolver<3>::generateR2()
{
	PetscErrorCode ierr;
	PetscInt       mstart, nstart, pstart, m, n, p, i, j, k;
	PetscReal      ***qx, ***qy, ***qz;
	PetscReal      ***bc2;
	Vec            bc2Global;

	ierr = VecSet(r2, 0.0); CHKERRV(ierr);
	ierr = DMCompositeGetAccess(lambdaPack, r2,  &bc2Global); CHKERRV(ierr);

	ierr = DMDAVecGetArray(uda, qxLocal, &qx); CHKERRV(ierr);
	ierr = DMDAVecGetArray(vda, qyLocal, &qy); CHKERRV(ierr);
	ierr = DMDAVecGetArray(wda, qzLocal, &qz); CHKERRV(ierr);

	ierr = DMDAVecGetArray(pda, bc2Global, &bc2); CHKERRV(ierr);
	ierr = DMDAGetCorners(pda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRV(ierr);
	// x-faces
	if(flowDesc->bc[0][XPLUS].type != PERIODIC) // don't update if the BC type is periodic
	{
		//loop over all points on the x-face
		for(k=pstart; k<pstart+p; k++)
		{
			for(j=nstart; j<nstart+n; j++)
			{
				// -X
				if(mstart == 0) // if the -X face is in the current process
				{
					bc2[k][j][0] -= qx[k][j][-1];
				}
				// +X
				if(mstart+m-1 == mesh->nx-1) // if the +X face is in the current process
				{
					bc2[k][j][mesh->nx-1] += qx[k][j][mesh->nx-1];
				}
			}
		}
	}
	// y-faces
	if(flowDesc->bc[1][YPLUS].type != PERIODIC) // don't update if the BC type is periodic
	{
		//loop over all points on the y-face
		for(k=pstart; k<pstart+p; k++)
		{
			for(i=mstart; i<mstart+m; i++)
			{
				// -Y
				if(nstart == 0) // if the -Y face is in the current process
				{
					bc2[k][0][i] -= qy[k][-1][i];
				}
				// +Y
				if(nstart+n-1 == mesh->ny-1) // if the +Y face is in the current process
				{
					bc2[k][mesh->ny-1][i] += qy[k][mesh->ny-1][i];
				}
			}
		}
	}
	// z-faces
	if(flowDesc->bc[2][ZPLUS].type != PERIODIC) // don't update if the BC type is periodic
	{
		//loop over all points on the z-face
		for(j=nstart; j<nstart+n; j++)
		{
			for(i=mstart; i<mstart+m; i++)
			{
				// -Z
				if(pstart == 0) // if the -Z face is in the current process
				{
					bc2[0][j][i] -= qz[-1][j][i];
				}
				// +Z
				if(pstart+p-1 == mesh->nz-1) // if the +Z face is in the current process
				{
					bc2[mesh->nz-1][j][i] += qz[mesh->nz-1][j][i];
				}
			}
		}
	}
	ierr = DMDAVecRestoreArray(pda, bc2Global, &bc2); CHKERRV(ierr);

	ierr = DMDAVecRestoreArray(uda, qxLocal, &qx); CHKERRV(ierr);
	ierr = DMDAVecRestoreArray(vda, qyLocal, &qy); CHKERRV(ierr);
	ierr = DMDAVecRestoreArray(wda, qzLocal, &qz); CHKERRV(ierr);

	ierr = DMCompositeRestoreAccess(lambdaPack, r2,  &bc2Global); CHKERRV(ierr);
}