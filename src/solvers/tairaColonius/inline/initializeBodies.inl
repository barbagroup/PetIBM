/*! Implementation of the method `initializeBodies` of the class `TairaColoniusSolver`.
 * \file initializeBodies.inl
 */


#include "yaml-cpp/yaml.h"


/*!
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


/*!
 * \brief Initializes the immersed boundaries.
 *
 * Parses the input file containing the list of the immersed boundaries using
 * YAML-CPP.
 */
// template <PetscInt dim>
// PetscErrorCode TairaColoniusSolver<dim>::initializeBodies()
// {
//   PetscErrorCode ierr;
//   PetscFunctionBeginUser;

//   std::string filePath = NavierStokesSolver<dim>::parameters->directory+"/bodies.yaml";
//   char path[PETSC_MAX_PATH_LEN];
//   PetscBool found;
//   PetscOptionsGetString(NULL, NULL, "-bodies", path, sizeof(path), &found);
//   if (found)
//     filePath = std::string(path);

//   YAML::Node nodes = YAML::LoadFile(filePath);
//   numBodies = nodes.size();
//   for (PetscInt i=0; i<numBodies; i++)
//   {
//     const YAML::Node &node = nodes[i];
//     std::string type = node["type"].as<std::string>();
//     if (type == "points")
//     {
//       size_t last = filePath.find_last_of("/");
//       std::string directory = filePath.substr(0, last+1);
//       std::string pointsFilePath = directory + node["pointsFile"].as<std::string>();
//       bodies.push_back(Body<dim>(pointsFilePath));
//     }
//   }

//   PetscFunctionReturn(0);
// } // initializeBodies
