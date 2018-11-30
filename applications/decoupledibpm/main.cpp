/**
 * \file decoupledibpm/main.cpp
 * \brief Main function of decoupled IBPM solver (Li et. al. 2016).
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 * \see decoupledibpm
 * \ingroup decoupledibpm
 */

#include <petscsys.h>
#include <yaml-cpp/yaml.h>

#include <petibm/parser.h>

#include "decoupledibpm.h"

/**
 * \defgroup decoupledibpm Decoupled IBPM solver (Li et. al. 2016)
 * \brief Implementation of parallel decoupled IBPM solver (Li et. al. 2016).
 *
 * This is an example of using PetIBM to build a parallel incompressible flow
 * solver with the immersed-boundary method proposed by Li et. al. 2016. We name
 * this solver \a decoupled \a IBPM".
 *
 * If readers are interested in using this solver instead of coding,
 * please refer to
 * \ref md_doc_markdowns_runpetibm "Running PetIBM",
 * \ref md_doc_markdowns_examples2d "2D Examples", and
 * \ref md_doc_markdowns_examples3d "3D Examples".
 *
 * \b Reference: \n
 * \li Li, R. Y., Xie, C. M., Huang, W. X., & Xu, C. X. (2016). An efficient
 * immersed boundary projection method for flow over complex/moving boundaries.
 * Computers & Fluids, 140, 122-135.
 *
 * \see nssolver, ibpm
 * \ingroup apps
 */

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    YAML::Node config;
    DecoupledIBPMSolver solver;

    ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);
    ierr = PetscLogDefaultBegin(); CHKERRQ(ierr);

    // parse configuration files; store info in YAML node
    ierr = petibm::parser::getSettings(config); CHKERRQ(ierr);

    // initialize the decoupled IBPM solver
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

    // destroy the decoupled IBPM solver
    ierr = solver.destroy(); CHKERRQ(ierr);

    ierr = PetscFinalize(); CHKERRQ(ierr);

    return 0;
}  // main
