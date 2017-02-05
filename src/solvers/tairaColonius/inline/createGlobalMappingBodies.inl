/*! Implementation of the method `createGlobalMappingBodies` of the class `TairaColoniusSolver`.
 * \file createGlobalMappingBodies.inl
 */


/*!
 * \brief Maps local to global indices for the Lagrangian forces.
 *
 * The mapping is stored in the Body objects.
 * Its maps the index of Lagrangian point to its global index in the vector lambda.
 */
template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::createGlobalMappingBodies()
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  // get the global offset index
  PetscInt lambdaStart;
  ierr = VecGetOwnershipRange(NavierStokesSolver<dim>::lambda, &lambdaStart, NULL); CHKERRQ(ierr);

  DMDALocalInfo pdaInfo, bdaInfo;
  ierr = DMDAGetLocalInfo(NavierStokesSolver<dim>::pda, &pdaInfo); CHKERRQ(ierr);

  PetscInt offset = lambdaStart + (pdaInfo.xm)*(pdaInfo.ym);
  for (auto &body : bodies)
  {
    // offset is passed by reference and updated in the Body method
    ierr = body.registerGlobalIdxPoints(offset); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
} // createGlobalMappingBodies
