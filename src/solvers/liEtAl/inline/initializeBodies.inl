/*! Implementation of the method `initializeBodies` of the class `LiEtAlSolver`.
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
PetscErrorCode LiEtAlSolver<dim>::initializeBodies()
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;

  std::string filePath = NavierStokesSolver<dim>::parameters->directory+"/bodies.yaml";
  // possibility to overwrite the path of the YAML configuration file
  // using the command-line parameter `-bodies <yaml-file-path>`
  char path[PETSC_MAX_PATH_LEN];
  PetscBool found;
  PetscOptionsGetString(NULL, NULL, "-bodies", path, sizeof(path), &found);
  if (found)
    filePath = std::string(path);

  YAML::Node nodes = YAML::LoadFile(filePath);
  numBodies = nodes.size();
  bodies.resize(numBodies);
  for (PetscInt i=0; i<numBodies; i++)
  {
    const YAML::Node &node = nodes[i];
    bodies[i] = Body<dim>();
    std::string type = node["type"].as<std::string>();
    if (type == "points")
    {
      std::string directory = NavierStokesSolver<dim>::parameters->directory;
      std::string pointsFileName = node["pointsFile"].as<std::string>();
      std::string pointsFilePath = directory + "/" + pointsFileName;
      ierr = bodies[i].readFromFile(pointsFilePath); CHKERRQ(ierr);
    }
    ierr = bodies[i].registerCellOwners(NavierStokesSolver<dim>::mesh); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
} // initializeBodies
