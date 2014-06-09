#include <petscdmcomposite.h>

template <>
PetscErrorCode NavierStokesSolver<2>::createVecs()
{
	PetscErrorCode    ierr;
	
	// local velocity fluxes
	ierr = DMCreateLocalVector(uda, &qxLocal); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(vda, &qyLocal); CHKERRQ(ierr);
	
	// global vectors
	ierr = DMCreateGlobalVector(qPack, &q); CHKERRQ(ierr); // velocity fluxes
	ierr = VecDuplicate(q, &qStar);        CHKERRQ(ierr); // intermediate velocity flux
	ierr = VecDuplicate(q, &H);            CHKERRQ(ierr); // convective term
	ierr = VecDuplicate(q, &rn);           CHKERRQ(ierr); // explicit terms
	ierr = VecDuplicate(q, &bc1);          CHKERRQ(ierr); // boundary conditions from implicit terms
	ierr = VecDuplicate(q, &rhs1);         CHKERRQ(ierr); // right-hand side for the intermediate-velocity solve
	ierr = VecDuplicate(q, &MHat);         CHKERRQ(ierr); // 
	ierr = VecDuplicate(q, &RInv);         CHKERRQ(ierr); // 
	ierr = VecDuplicate(q, &BN);           CHKERRQ(ierr); // 
	ierr = VecDuplicate(q, &temp);         CHKERRQ(ierr); // 

	ierr = DMCreateGlobalVector(lambdaPack, &lambda); CHKERRQ(ierr); // pressure
	ierr = VecDuplicate(lambda, &r2);              CHKERRQ(ierr); // 
	ierr = VecDuplicate(lambda, &rhs2);            CHKERRQ(ierr); // 

	return 0;
}

template <>
PetscErrorCode NavierStokesSolver<3>::createVecs()
{
	PetscErrorCode    ierr;
	
	// local velocity fluxes
	ierr = DMCreateLocalVector(uda, &qxLocal); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(vda, &qyLocal); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(wda, &qzLocal); CHKERRQ(ierr);
	
	// global vectors
	ierr = DMCreateGlobalVector(qPack, &q); CHKERRQ(ierr); // velocity fluxes
	ierr = VecDuplicate(q, &qStar);        CHKERRQ(ierr); // intermediate velocity flux
	ierr = VecDuplicate(q, &H);            CHKERRQ(ierr); // convective term
	ierr = VecDuplicate(q, &rn);           CHKERRQ(ierr); // explicit terms
	ierr = VecDuplicate(q, &bc1);          CHKERRQ(ierr); // boundary conditions from implicit terms
	ierr = VecDuplicate(q, &rhs1);         CHKERRQ(ierr); // right-hand side for the intermediate-velocity solve
	ierr = VecDuplicate(q, &MHat);         CHKERRQ(ierr); // 
	ierr = VecDuplicate(q, &RInv);         CHKERRQ(ierr); // 
	ierr = VecDuplicate(q, &BN);           CHKERRQ(ierr); // 
	ierr = VecDuplicate(q, &temp);         CHKERRQ(ierr); // 

	ierr = DMCreateGlobalVector(lambdaPack, &lambda); CHKERRQ(ierr); // pressure
	ierr = VecDuplicate(lambda, &r2);              CHKERRQ(ierr); // 
	ierr = VecDuplicate(lambda, &rhs2);            CHKERRQ(ierr); // 

	return 0;
}
