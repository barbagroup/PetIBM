/***************************************************************************//**
 * \file setNullSpace.inl
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the method `setNullSpace` of the class `NavierStokesSolver`.
 */


/**
 * \brief Sets the nullspace for the Poisson system.
 * 
 * The null space for the Poisson system in a flow with no immersed boundary
 * is the constant vector of size equal to the number of pressure variables. 
 * Such a null space can be specified and automatically handled by PETSc by 
 * passing PETSC_TRUE as the second parameter in the function MatNullSpaceCreate.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::setNullSpace()
{
	PetscErrorCode ierr;
	
	MatNullSpace nsp;
	ierr = MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, NULL, &nsp); CHKERRQ(ierr);
	ierr = KSPSetNullSpace(ksp2, nsp); CHKERRQ(ierr);
	ierr = MatNullSpaceDestroy(&nsp); CHKERRQ(ierr);

	return 0;
} // setNullSpace