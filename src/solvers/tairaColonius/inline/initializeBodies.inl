/*! Implementation of the method `initializeBodies` of the class `TairaColoniusSolver`.
 * \file initializeBodies.inl
 */


#include "yaml-cpp/yaml.h"


/*!
 * \brief Initializes the immersed boundaries.
 *
 * Parses the input file containing the list of the immersed boundaries using
 * YAML-CPP.
 */
template <PetscInt dim>
PetscErrorCode TairaColoniusSolver<dim>::initializeBodies()
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;

  std::string filePath = NavierStokesSolver<dim>::parameters->directory+"/bodies.yaml";
  char path[PETSC_MAX_PATH_LEN];
  PetscBool found;
  PetscOptionsGetString(NULL, NULL, "-bodies", path, sizeof(path), &found);
  if (found)
    filePath = std::string(path);

  YAML::Node nodes = YAML::LoadFile(filePath);
  numBodies = nodes.size();
  numLagrangianPoints = 0;
  bodies.resize(numBodies);
  for (PetscInt i=0; i<numBodies; i++)
  {
    const YAML::Node &node = nodes[i];
    bodies[i] = Body<dim>();
    ierr = bodies[i].initialize(); CHKERRQ(ierr);
    std::string type = node["type"].as<std::string>();
    if (type == "points")
    {
      size_t last = filePath.find_last_of("/");
      std::string directory = filePath.substr(0, last+1);
      std::string pointsFilePath = directory + node["pointsFile"].as<std::string>();
      ierr = bodies[i].readFromFile(pointsFilePath); CHKERRQ(ierr);
    }
    numLagrangianPoints += bodies[i].numPoints;
    ierr = bodies[i].setCellOwners(NavierStokesSolver<dim>::mesh); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
} // initializeBodies
