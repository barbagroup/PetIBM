/***************************************************************************//**
* Allocate memory for the vectors required in the simulation using the 
* corresponding distributed array structures.
*/
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::createVecs()
{
	PetscErrorCode    ierr;
	
	// local vectors to store velocity fluxes
	ierr = DMCreateLocalVector(uda, &qxLocal); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(vda, &qyLocal); CHKERRQ(ierr);
	if(dim==3)
	{
		ierr = DMCreateLocalVector(wda, &qzLocal); CHKERRQ(ierr);
	}
	
	// global vectors
	ierr = DMCreateGlobalVector(qPack, &q); CHKERRQ(ierr); // velocity fluxes
	ierr = VecDuplicate(q, &qStar);        CHKERRQ(ierr); // intermediate velocity flux
	ierr = VecDuplicate(q, &H);            CHKERRQ(ierr); // convective term
	ierr = VecDuplicate(q, &rn);           CHKERRQ(ierr); // explicit terms
	ierr = VecDuplicate(q, &bc1);          CHKERRQ(ierr); // boundary conditions from implicit terms
	ierr = VecDuplicate(q, &rhs1);         CHKERRQ(ierr); // right-hand side for the intermediate-velocity solve
	ierr = VecDuplicate(q, &MHat);         CHKERRQ(ierr); // 
	ierr = VecDuplicate(q, &RInv);         CHKERRQ(ierr); // 
	ierr = VecDuplicate(q, &BN);           CHKERRQ(ierr); // approximate inverse of `A`
	ierr = VecDuplicate(q, &temp);         CHKERRQ(ierr); // 

	ierr = DMCreateGlobalVector(lambdaPack, &lambda); CHKERRQ(ierr); // pressure
	ierr = VecDuplicate(lambda, &r2);              CHKERRQ(ierr); // 
	ierr = VecDuplicate(lambda, &rhs2);            CHKERRQ(ierr); // right-hand size for the Poisson solve

	return 0;
}