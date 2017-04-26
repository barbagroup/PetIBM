/*! Implementation of the method `createGlobalMappingBodies` of the class `LiEtAlSolver`.
 * \file createGlobalMappingBodies.inl
 */


/*!
 * \brief Maps local to global indices for the Lagrangian forces.
 *
 * The mapping is stored in the Body objects.
 * Its maps the index of Lagrangian point to its global index in the vector f.
 */
template <PetscInt dim>
PetscErrorCode LiEtAlSolver<dim>::createGlobalMappingBodies()
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  // get the global offset index
  PetscInt fTildeStart;
  ierr = VecGetOwnershipRange(fTilde, &fTildeStart, NULL); CHKERRQ(ierr);

  PetscInt offset = fTildeStart;
  for (auto &body : bodies)
  {
    // offset is passed by reference and updated in the Body method
    ierr = body.registerGlobalIdxPoints(offset); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
} // createGlobalMappingBodies
