template <>
void NavierStokesSolver<2>::generateM()
{
	PetscErrorCode ierr;
	PetscInt       mstart, nstart, m, n, i, j;
	Vec            MxGlobal, MyGlobal;
	PetscReal      **Mx, **My;

	ierr = DMCompositeGetAccess(pack, M, &MxGlobal, &MyGlobal); CHKERRV(ierr);

	// x-direction
	ierr = DMDAVecGetArray(uda, MxGlobal, &Mx); CHKERRV(ierr);
	ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRV(ierr);
	for(j=nstart; j<nstart+n; j++)
	{
		for(i=mstart; i<mstart+m; i++)
		{
			Mx[j][i] = 0.5*(mesh->dx[i] + mesh->dx[i+1]);
		}
	}
	ierr = DMDAVecRestoreArray(uda, MxGlobal, &Mx); CHKERRV(ierr);

	// y-direction
	ierr = DMDAVecGetArray(vda, MyGlobal, &My); CHKERRV(ierr);
	ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRV(ierr);
	for(j=nstart; j<nstart+n; j++)
	{
		for(i=mstart; i<mstart+m; i++)
		{
			My[j][i] = 0.5*(mesh->dy[j] + mesh->dy[j+1]);
		}
	}
	ierr = DMDAVecRestoreArray(vda, MyGlobal, &My); CHKERRV(ierr);

	ierr = DMCompositeRestoreAccess(pack, M, &MxGlobal, &MyGlobal); CHKERRV(ierr);
}

template <>
void NavierStokesSolver<3>::generateM()
{

}