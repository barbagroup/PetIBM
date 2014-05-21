#include <petscdmcomposite.h>

template <>
void NavierStokesSolver<2>::createVecs()
{
	PetscErrorCode    ierr;
	
	// local velocity fluxes
	ierr = DMCreateLocalVector(uda, &qxLocal); CHKERRV(ierr);
	ierr = DMCreateLocalVector(vda, &qyLocal); CHKERRV(ierr);
	
	// global vectors
	ierr = DMCreateGlobalVector(qPack, &q); CHKERRV(ierr); // velocity fluxes
	ierr = VecDuplicate(q, &qStar);        CHKERRV(ierr); // intermediate velocity flux
	ierr = VecDuplicate(q, &H);            CHKERRV(ierr); // convective term
	ierr = VecDuplicate(q, &rn);           CHKERRV(ierr); // explicit terms
	ierr = VecDuplicate(q, &bc1);          CHKERRV(ierr); // boundary conditions from implicit terms
	ierr = VecDuplicate(q, &rhs1);         CHKERRV(ierr); // right-hand side for the intermediate-velocity solve
	ierr = VecDuplicate(q, &MHat);         CHKERRV(ierr); // 
	ierr = VecDuplicate(q, &RInv);         CHKERRV(ierr); // 
	ierr = VecDuplicate(q, &BN);           CHKERRV(ierr); // 
	ierr = VecDuplicate(q, &temp);         CHKERRV(ierr); // 

	ierr = DMCreateGlobalVector(phiPack, &phi); CHKERRV(ierr); // pressure
	ierr = VecDuplicate(phi, &r2);              CHKERRV(ierr); // 
	ierr = VecDuplicate(phi, &rhs2);            CHKERRV(ierr); // 
}

template <>
void NavierStokesSolver<3>::createVecs()
{
	PetscErrorCode    ierr;
	
	// local velocity fluxes
	ierr = DMCreateLocalVector(uda, &qxLocal); CHKERRV(ierr);
	ierr = DMCreateLocalVector(vda, &qyLocal); CHKERRV(ierr);
	ierr = DMCreateLocalVector(wda, &qzLocal); CHKERRV(ierr);
	
	// global vectors
	ierr = DMCreateGlobalVector(qPack, &q); CHKERRV(ierr); // velocity fluxes
	ierr = VecDuplicate(q, &qStar);        CHKERRV(ierr); // intermediate velocity flux
	ierr = VecDuplicate(q, &H);            CHKERRV(ierr); // convective term
	ierr = VecDuplicate(q, &rn);           CHKERRV(ierr); // explicit terms
	ierr = VecDuplicate(q, &bc1);          CHKERRV(ierr); // boundary conditions from implicit terms
	ierr = VecDuplicate(q, &rhs1);         CHKERRV(ierr); // right-hand side for the intermediate-velocity solve
	ierr = VecDuplicate(q, &MHat);         CHKERRV(ierr); // 
	ierr = VecDuplicate(q, &RInv);         CHKERRV(ierr); // 
	ierr = VecDuplicate(q, &BN);           CHKERRV(ierr); // 
	ierr = VecDuplicate(q, &temp);         CHKERRV(ierr); // 

	ierr = DMCreateGlobalVector(phiPack, &phi); CHKERRV(ierr); // pressure
	ierr = VecDuplicate(phi, &r2);              CHKERRV(ierr); // 
	ierr = VecDuplicate(phi, &rhs2);            CHKERRV(ierr); // 
}
