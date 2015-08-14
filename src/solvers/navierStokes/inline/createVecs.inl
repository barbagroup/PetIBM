/***************************************************************************//**
 * \file createVecs.inl
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the method `createVecs` of the class `NavierStokesSolver`.
 */


/**
 * \brief Allocates memory for the vectors required in the simulation using the 
 *        corresponding distributed array structures.
 */
template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::createVecs()
{
  PetscErrorCode ierr;
  
  // create local vectors to store fluxes in each direction
  ierr = DMCreateLocalVector(uda, &qxLocal); CHKERRQ(ierr);
  ierr = DMCreateLocalVector(vda, &qyLocal); CHKERRQ(ierr);
  if (dim == 3)
  {  
    ierr = DMCreateLocalVector(wda, &qzLocal); CHKERRQ(ierr);
  }
  // create global vectors needed for velocity system
  ierr = DMCreateGlobalVector(qPack, &q); CHKERRQ(ierr); // velocity fluxes
  ierr = VecDuplicate(q, &qStar); CHKERRQ(ierr);         // intermediate velocity fluxes
  ierr = VecDuplicate(q, &H); CHKERRQ(ierr);             // convective terms
  ierr = VecDuplicate(q, &rn); CHKERRQ(ierr);            // explicit terms
  ierr = VecDuplicate(q, &bc1); CHKERRQ(ierr);           // boundary conditions from implicit terms
  ierr = VecDuplicate(q, &rhs1); CHKERRQ(ierr);          // RHS of intermediate-velocity system
  ierr = VecDuplicate(q, &MHat); CHKERRQ(ierr);          // scaling mass matrix 
  ierr = VecDuplicate(q, &RInv); CHKERRQ(ierr);          // scaling matrix 
  ierr = VecDuplicate(q, &BN); CHKERRQ(ierr);            // approximate inverse of `A`
  ierr = VecDuplicate(q, &temp); CHKERRQ(ierr);          // temporary vector 
  // create global vectors needed for Poisson system
  ierr = DMCreateGlobalVector(lambdaPack, &lambda); CHKERRQ(ierr); // pressure
  ierr = VecDuplicate(lambda, &r2); CHKERRQ(ierr);                 //  
  ierr = VecDuplicate(lambda, &rhs2); CHKERRQ(ierr);               // right-hand size for the Poisson solve

  return 0;
} // createVecs