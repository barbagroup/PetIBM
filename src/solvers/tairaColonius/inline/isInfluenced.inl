/***************************************************************************//**
 * \file isInfluenced.inl
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the method `isInfluenced` of the class `TairaColoniusSolver`.
 */


/**
 * \brief Determines of an Eulerian grid point is in the disk of influence of a
 *        Lagrangian body point.
 */
template <PetscInt dim>
PetscBool TairaColoniusSolver<dim>::isInfluenced(PetscReal xGrid, PetscReal yGrid, PetscReal xBody, PetscReal yBody, PetscReal radius, PetscReal *disp)
{
	PetscReal width[2];
	PetscReal nx = NavierStokesSolver<dim>::mesh->nx,
	          ny = NavierStokesSolver<dim>::mesh->ny;
	
	std::vector<PetscReal> &xMesh = NavierStokesSolver<dim>::mesh->x,
	                       &yMesh = NavierStokesSolver<dim>::mesh->y;
	
	width[0] = xMesh[nx]-xMesh[0];
	width[1] = yMesh[ny]-yMesh[0];

	disp[0] = fabs(xGrid - xBody);
	disp[1] = fabs(yGrid - yBody);

	if (NavierStokesSolver<dim>::flow->boundaries[XPLUS][0].type == PERIODIC && disp[0] > width[0]-disp[0])
		disp[0] = width[0] - disp[0];
	if (NavierStokesSolver<dim>::flow->boundaries[YPLUS][0].type == PERIODIC && disp[1] > width[1]-disp[1])
		disp[1] = width[1] - disp[1];

	return (disp[0] < radius && disp[1] < radius)? PETSC_TRUE : PETSC_FALSE;
} // isInfluenced


/**
 * \brief Determines of an Eulerian grid point is in the sphere of influence of a
 *        Lagrangian body point.
 */
template <PetscInt dim>
PetscBool TairaColoniusSolver<dim>::isInfluenced(PetscReal xGrid, PetscReal yGrid, PetscReal zGrid, PetscReal xBody, PetscReal yBody, PetscReal zBody, PetscReal radius, PetscReal *disp)
{
	PetscReal width[3];
	PetscReal nx = NavierStokesSolver<dim>::mesh->nx,
	          ny = NavierStokesSolver<dim>::mesh->ny,
	          nz = NavierStokesSolver<dim>::mesh->nz;

	std::vector<PetscReal> &xMesh = NavierStokesSolver<dim>::mesh->x,
	                       &yMesh = NavierStokesSolver<dim>::mesh->y,
	                       &zMesh = NavierStokesSolver<dim>::mesh->z;
	
	width[0] = xMesh[nx]-xMesh[0];
	width[1] = yMesh[ny]-yMesh[0];
	width[2] = zMesh[nz]-zMesh[0];

	disp[0] = fabs(xGrid - xBody);
	disp[1] = fabs(yGrid - yBody);
	disp[2] = fabs(zGrid - zBody);

	if (NavierStokesSolver<dim>::flow->boundaries[XPLUS][0].type == PERIODIC && disp[0] > width[0]-disp[0])
		disp[0] = width[0] - disp[0];
	if (NavierStokesSolver<dim>::flow->boundaries[YPLUS][0].type == PERIODIC && disp[1] > width[1]-disp[1])
		disp[1] = width[1] - disp[1];
	if (NavierStokesSolver<dim>::flow->boundaries[ZPLUS][0].type == PERIODIC && disp[2] > width[2]-disp[2])
		disp[2] = width[2] - disp[2];

	return (disp[0] < radius && disp[1] < radius && disp[2] < radius)? PETSC_TRUE : PETSC_FALSE;
} // isInfluenced