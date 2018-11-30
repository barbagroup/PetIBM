/**
 * \file ibpm/main.cpp
 * \brief Main function of IBPM solver (Taira & Colonius 2007).
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 * \see ibpm
 * \ingroup ibpm
 */

#include <petscsys.h>
#include <yaml-cpp/yaml.h>

#include <petibm/parser.h>

#include "ibpm.h"

/**
 * \defgroup ibpm IBPM solver (Taira & Colonius 2007)
 * \brief Implementation of parallel IBPM solver (Taira & Colonius 2007).
 *
 * This is an example of using PetIBM to build a parallel incompressible flow
 * solver with the immersed-boundary method proposed by Taira & Colonius 2007.
 * We name this solver \a IBPM.
 *
 * If readers are interested in using this solver instead of coding,
 * please refer to
 * \ref md_doc_markdowns_runpetibm "Running PetIBM",
 * \ref md_doc_markdowns_examples2d "2D Examples", and
 * \ref md_doc_markdowns_examples3d "3D Examples".
 *
 * \b Reference: \n
 * \li Taira, K., & Colonius, T. (2007). The immersed boundary method: a
 * projection approach. Journal of Computational Physics, 225(2), 2118-2137.
 *
 * \see nssolver, decoupledibpm
 * \ingroup apps
 */

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    YAML::Node config;
    IBPMSolver solver;

    ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);
    ierr = PetscLogDefaultBegin(); CHKERRQ(ierr);

    // parse configuration files; store info in YAML node
    ierr = petibm::parser::getSettings(config); CHKERRQ(ierr);

    // initialize the IBPM solver
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

    // destroy the IBPM solver
    ierr = solver.destroy(); CHKERRQ(ierr);

    ierr = PetscFinalize(); CHKERRQ(ierr);

    return 0;
}  // main
