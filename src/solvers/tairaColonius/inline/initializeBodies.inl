/***************************************************************************//**
 * \file initializeBodies.inl
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the method `initializeBodies` of the class `TairaColoniusSolver`.
 */


/**
 * \brief Initializes the immersed boundaries.
 */
template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::initializeBodies()
{
  numBodies = 1;
  bodies.resize(numBodies);
  for (PetscInt l=0; l<numBodies; l++)
  {
    bodies[l] = Body<dim>(NavierStokesSolver<dim>::parameters->directory+"/bodies.yaml");
  }

  return 0;
} // initializeBodies