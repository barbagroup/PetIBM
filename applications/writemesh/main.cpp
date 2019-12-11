/**
 * \file mesh/main.cpp
 * \brief Small application to create and write mesh coordinates to a file.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 * \see writemesh
 * \ingroup writemesh
 */

#include <petscsys.h>
#include <yaml-cpp/yaml.h>

#include <petibm/mesh.h>
#include <petibm/parser.h>

/**
 * \defgroup writemesh Pre-processing utility: writemesh
 * \brief A pre-processing utility that creates a mesh and writes
 *        the gridline coordinates to a HDF5 file.
 *
 * \ingroup apps
 */

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    YAML::Node config;
    petibm::type::Mesh mesh;

    ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);
    ierr = PetscLogDefaultBegin(); CHKERRQ(ierr);

    // parse configuration files; store info in YAML node
    ierr = petibm::parser::getSettings(config); CHKERRQ(ierr);

    // create mesh
    ierr = petibm::mesh::createMesh(PETSC_COMM_WORLD, config, mesh);
    CHKERRQ(ierr);

    // write grid data into HDF5 file
    std::string filePath = config["output"].as<std::string>() + "/grid.h5";
    char s[PETSC_MAX_PATH_LEN];
    PetscBool flag = PETSC_FALSE;
    ierr = PetscOptionsGetString(nullptr, nullptr, "-file", s, sizeof(s),
                                 &flag); CHKERRQ(ierr);
    if (flag) filePath = s;
    ierr = mesh->write(filePath); CHKERRQ(ierr);

    ierr = mesh->destroy(); CHKERRQ(ierr);
    ierr = PetscFinalize(); CHKERRQ(ierr);

    return 0;
}  // main
