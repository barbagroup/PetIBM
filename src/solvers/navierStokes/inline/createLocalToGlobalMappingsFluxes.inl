/***************************************************************************//**
 * \file createLocalToGlobalMappingsFluxes.inl
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the method `createLocalToGlobalMappingsFluxes` 
 *        of the class `NavierStokesSolver`.
 */


/**
 * \brief Maps local multi-dimensional indices to global index for the fluxes.
 * 
 * Vectors stored as distributed arrays can be accessed using multi-dimensional 
 * arrays on every process, with each index referring to the numbering along 
 * each cartesian direction. The elements of the vector also have a global
 * ordering. This function generates the map from the multi-dimensional
 * indexing of each of the local flux vectors `qx`, `qy` and `qz`, to the global 
 * indices of the composite flux vector `q`.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::createLocalToGlobalMappingsFluxes()
{
	return 0;
} // createLocalToGlobalMappingsFluxes


// two-dimensional specialization
template <>
PetscErrorCode NavierStokesSolver<2>::createLocalToGlobalMappingsFluxes()
{
	PetscErrorCode ierr;
	PetscInt i, j,           // loop indices
					 m, n,           // local number of nodes along each direction
					 mstart, nstart; // starting indices

	// get global index of first local element of fluxes
	PetscInt globalIdx;
	ierr = VecGetOwnershipRange(q, &globalIdx, NULL); CHKERRQ(ierr);

	// populate local vectors with the global indices
	// set value to -1 if the cell is outside the domain
	// fluxes in x-direction
	ierr = DMCreateLocalVector(uda, &uMapping); CHKERRQ(ierr);
	PetscReal **uMappingArray;
	ierr = DMDAVecGetArray(uda, uMapping, &uMappingArray); CHKERRQ(ierr);
	ierr = DMDAGetGhostCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
	for (j=nstart; j<nstart+n; j++)
	{
		for (i=mstart; i<mstart+m; i++)
		{
			uMappingArray[j][i] = -1;
			if (i > mstart && i < mstart+m-1 && j > nstart && j < nstart+n-1)
			{
				uMappingArray[j][i] = globalIdx;
				globalIdx++;
			}
		}
	}
	ierr = DMDAVecRestoreArray(uda, uMapping, &uMappingArray); CHKERRQ(ierr);
	// fluxes in y-direction
	ierr = DMCreateLocalVector(vda, &vMapping); CHKERRQ(ierr);
	PetscReal **vMappingArray;
	ierr = DMDAVecGetArray(vda, vMapping, &vMappingArray); CHKERRQ(ierr);
	ierr = DMDAGetGhostCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
	for (j=nstart; j<nstart+n; j++)
	{
		for (i=mstart; i<mstart+m; i++)
		{
			vMappingArray[j][i] = -1;
			if (i > mstart && i < mstart+m-1 && j > nstart && j < nstart+n-1)
			{
				vMappingArray[j][i] = globalIdx;
				globalIdx++;
			}	
		}
	}
	ierr = DMDAVecRestoreArray(vda, vMapping, &vMappingArray); CHKERRQ(ierr);

	// scatter from local to local to obtain correct values in ghost cells
	// fluxes in x-direction
	ierr = DMLocalToLocalBegin(uda, uMapping, INSERT_VALUES, uMapping); CHKERRQ(ierr);
	ierr = DMLocalToLocalEnd(uda, uMapping, INSERT_VALUES, uMapping); CHKERRQ(ierr);
	// fluxes in y-direction
	ierr = DMLocalToLocalBegin(vda, vMapping, INSERT_VALUES, vMapping); CHKERRQ(ierr);
	ierr = DMLocalToLocalEnd(vda, vMapping, INSERT_VALUES, vMapping); CHKERRQ(ierr);

	return 0;
} // createLocalToGlobalMappingsFluxes


// three-dimensional specialization
template <>
PetscErrorCode NavierStokesSolver<3>::createLocalToGlobalMappingsFluxes()		
{
	PetscErrorCode ierr;
	PetscInt i, j, k,                // loop indices
					 m, n, p,                // local number of nodes along each direction
					 mstart, nstart, pstart; // starting indices

	// get global index of first local element of fluxes
	PetscInt globalIdx;
	ierr = VecGetOwnershipRange(q, &globalIdx, NULL); CHKERRQ(ierr);

	// populate local vectors with the global indices
	// set value to -1 if the cell is outside the domain
	// fluxes in x-direction
	ierr = DMCreateLocalVector(uda, &uMapping); CHKERRQ(ierr);
	PetscReal ***uMappingArray;
	ierr = DMDAVecGetArray(uda, uMapping, &uMappingArray); CHKERRQ(ierr);
	ierr = DMDAGetGhostCorners(uda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
	for (k=pstart; k<pstart+p; k++)
	{
		for (j=nstart; j<nstart+n; j++)
		{
			for (i=mstart; i<mstart+m; i++)
			{
				uMappingArray[k][j][i] = -1;
				if (i > mstart && i < mstart+m-1 && j > nstart && j < nstart+n-1 && k > pstart && k < pstart+p-1)
				{
					uMappingArray[k][j][i] = globalIdx;
					globalIdx++;
				}
			}
		}
	}
	ierr = DMDAVecRestoreArray(uda, uMapping, &uMappingArray); CHKERRQ(ierr);
	// fluxes in y-direction
	ierr = DMCreateLocalVector(vda, &vMapping); CHKERRQ(ierr);
	PetscReal ***vMappingArray;
	ierr = DMDAVecGetArray(vda, vMapping, &vMappingArray); CHKERRQ(ierr);
	ierr = DMDAGetGhostCorners(vda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
	for (k=pstart; k<pstart+p; k++)
	{
		for (j=nstart; j<nstart+n; j++)
		{
			for (i=mstart; i<mstart+m; i++)
			{
				vMappingArray[k][j][i] = -1;
				if (i > mstart && i < mstart+m-1 && j > nstart && j < nstart+n-1 && k > pstart && k < pstart+p-1)
				{
					vMappingArray[k][j][i] = globalIdx;
					globalIdx++;
				}
			}
		}
	}
	ierr = DMDAVecRestoreArray(vda, vMapping, &vMappingArray); CHKERRQ(ierr);
	// fluxes in z-direction
	ierr = DMCreateLocalVector(wda, &wMapping); CHKERRQ(ierr);
	PetscReal ***wMappingArray;
	ierr = DMDAVecGetArray(wda, wMapping, &wMappingArray); CHKERRQ(ierr);
	ierr = DMDAGetGhostCorners(wda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
	for (k=pstart; k<pstart+p; k++)
	{
		for (j=nstart; j<nstart+n; j++)
		{
			for (i=mstart; i<mstart+m; i++)
			{
				wMappingArray[k][j][i] = -1;
				if (i > mstart && i < mstart+m-1 && j > nstart && j < nstart+n-1 && k > pstart && k < pstart+p-1)
				{
					wMappingArray[k][j][i] = globalIdx;
					globalIdx++;
				}
			}
		}
	}
	ierr = DMDAVecRestoreArray(wda, wMapping, &wMappingArray); CHKERRQ(ierr);

	// scatter from local to local to obtain correct values in ghost cells
	// fluxes in x-direction
	ierr = DMLocalToLocalBegin(uda, uMapping, INSERT_VALUES, uMapping); CHKERRQ(ierr);
	ierr = DMLocalToLocalEnd(uda, uMapping, INSERT_VALUES, uMapping); CHKERRQ(ierr);
	// fluxes in y-direction
	ierr = DMLocalToLocalBegin(vda, vMapping, INSERT_VALUES, vMapping); CHKERRQ(ierr);
	ierr = DMLocalToLocalEnd(vda, vMapping, INSERT_VALUES, vMapping); CHKERRQ(ierr);
	// fluxes in z-direction
	ierr = DMLocalToLocalBegin(wda, wMapping, INSERT_VALUES, wMapping); CHKERRQ(ierr);
	ierr = DMLocalToLocalEnd(wda, wMapping, INSERT_VALUES, wMapping); CHKERRQ(ierr);

	return 0;
} // createLocalToGlobalMappingsFluxes