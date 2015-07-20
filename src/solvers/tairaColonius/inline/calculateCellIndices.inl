template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::calculateCellIndices()
{
	std::vector<PetscReal> &xMesh = NavierStokesSolver<dim>::mesh->x,
	                       &yMesh = NavierStokesSolver<dim>::mesh->y,
	                       &zMesh = NavierStokesSolver<dim>::mesh->z;

	I.reserve(x.size());
	J.reserve(x.size());
	if(dim==3) K.reserve(x.size());

	PetscInt i=0, j=0, k=0;

	/// find the cell for the zeroth point
	while(xMesh[i+1] < x[0]) i++;
	while(yMesh[j+1] < y[0]) j++;
	if(dim==3)
	{
		while(zMesh[k+1] < z[0]) k++;
	}

	I.push_back(i);
	J.push_back(j);
	if(dim==3) K.push_back(k);

	for(size_t l=1; l<x.size(); l++)
	{
		// if the next boundary point is to the left of the current boundary point
		if(x[l] < x[l-1])
		{
			while(xMesh[i] > x[l])
				i--;
		}
		// if the next boundary point is to the right of the current boundary point
		else
		{
			while(xMesh[i+1] < x[l])
				i++;
		}
		// if the next boundary point is below the current boundary point
		if(y[l] < y[l-1])
		{
			while(yMesh[j] > y[l])
				j--;
		}
		// if the next boundary point is above the current boundary point
		else
		{
			while(yMesh[j+1] < y[l])
				j++;
		}
		if(dim==3)
		{
			// if the next boundary point is below the current boundary point
			if(z[l] < z[l-1])
			{
				while(zMesh[k] > z[l])
					k--;
			}
			// if the next boundary point is above the current boundary point
			else
			{
				while(zMesh[k+1] < z[l])
					k++;
			}
		}

		I.push_back(i);
		J.push_back(j);
		if(dim==3) K.push_back(k);
	}

	return 0;
} // calculateCellIndices