#include <petscdmda.h>
template <>
void NavierStokesSolver<2>::createLocalToGlobalMappings()
{
	PetscErrorCode ierr;
	PetscInt       m, n, i, j, mstart, nstart;
	Vec            globalX, globalY;
	PetscReal      **gx, **gy, **lx, **ly;
	PetscInt       localIdx;

	ierr = DMCreateGlobalVector(pack, &globalIndices); CHKERRV(ierr);
	ierr = DMCreateLocalVector(uda, &uMapping); CHKERRV(ierr);
	ierr = DMCreateLocalVector(vda, &vMapping); CHKERRV(ierr);

	// initialise global vector
	ierr = VecGetOwnershipRange(globalIndices, &localIdx, NULL); CHKERRV(ierr);
	ierr = DMCompositeGetAccess(pack, globalIndices, &globalX, &globalY); CHKERRV(ierr);

	// U
	ierr = DMDAVecGetArray(uda, globalX, &gx); CHKERRV(ierr);
	ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRV(ierr);
	for(j=nstart; j<nstart+n; j++)
	{
		for(i=mstart; i<mstart+m; i++)
		{
			gx[j][i] = localIdx;
			localIdx++;
		}
	}
	ierr = DMDAVecRestoreArray(uda, globalX, &gx); CHKERRV(ierr);

	// V
	ierr = DMDAVecGetArray(vda, globalY, &gy); CHKERRV(ierr);
	ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRV(ierr);
	for(j=nstart; j<nstart+n; j++)
	{
		for(i=mstart; i<mstart+m; i++)
		{
			gy[j][i] = localIdx;
			localIdx++;
		}
	}
	ierr = DMDAVecRestoreArray(vda, globalY, &gy); CHKERRV(ierr);

	// initialise local vector
	// set values to -1
	// U
	ierr = DMDAVecGetArray(uda, uMapping, &lx); CHKERRV(ierr);
	ierr = DMDAGetGhostCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRV(ierr);
	for(j=nstart; j<nstart+n; j++)
	{
		for(i=mstart; i<mstart+m; i++)
		{
			lx[j][i] = -1;
		}
	}
	ierr = DMDAVecRestoreArray(uda, uMapping, &lx); CHKERRV(ierr);

	// V
	ierr = DMDAVecGetArray(vda, vMapping, &ly); CHKERRV(ierr);
	ierr = DMDAGetGhostCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRV(ierr);
	for(j=nstart; j<nstart+n; j++)
	{
		for(i=mstart; i<mstart+m; i++)
		{
			ly[j][i] = -1;
		}
	}
	ierr = DMDAVecRestoreArray(vda, vMapping, &ly); CHKERRV(ierr);

	// scatter from global to local vectors
	ierr = DMGlobalToLocalBegin(uda, globalX, INSERT_VALUES, uMapping); CHKERRV(ierr);
	ierr = DMGlobalToLocalEnd(uda, globalX, INSERT_VALUES, uMapping); CHKERRV(ierr);

	ierr = DMGlobalToLocalBegin(vda, globalY, INSERT_VALUES, vMapping); CHKERRV(ierr);
	ierr = DMGlobalToLocalEnd(vda, globalY, INSERT_VALUES, vMapping); CHKERRV(ierr);

	ierr = DMCompositeRestoreAccess(pack, globalIndices, &globalX, &globalY); CHKERRV(ierr);

	ierr = VecDestroy(&globalIndices);
}

template <>
void NavierStokesSolver<3>::createLocalToGlobalMappings()		
{

}