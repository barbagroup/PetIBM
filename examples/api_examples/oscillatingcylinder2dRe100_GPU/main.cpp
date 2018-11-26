/**
 * \file main.cpp
 * \brief Main function for the simulation of a two-dimensional
 *        inline-oscillating cylinder.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 */

#include <petscsys.h>
#include <yaml-cpp/yaml.h>

#include <petibm/parser.h>

#include "oscillatingcylinder.h"

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    YAML::Node config;
    OscillatingCylinderSolver solver;

    ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);
    ierr = PetscLogDefaultBegin(); CHKERRQ(ierr);

    // parse configuration files; store info in YAML node
    ierr = petibm::parser::getSettings(config); CHKERRQ(ierr);

    // initialize the solver
    ierr = solver.init(PETSC_COMM_WORLD, config); CHKERRQ(ierr);
    ierr = solver.ioInitialData(); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,
                       "Completed initialization stage\n"); CHKERRQ(ierr);

    // integrate the solution in time
    while (!solver.finished())
    {
        // compute the solution at the next time step
        ierr = solver.advance(); CHKERRQ(ierr);
        // output data to files
        ierr = solver.write(); CHKERRQ(ierr);
    }

    // destroy the solver
    ierr = solver.destroy(); CHKERRQ(ierr);

    ierr = PetscFinalize(); CHKERRQ(ierr);

    return 0;
}  // main
