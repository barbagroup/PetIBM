/**
 * \file navierstokes/main.cpp
 * \brief Main function of the Navier-Stokes solver.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 * \see nssolver
 * \ingroup nssolver
 */

#include <petscsys.h>
#include <yaml-cpp/yaml.h>

#include <petibm/parser.h>

#include "navierstokes.h"

/**
 * \defgroup apps Flow solvers and utilities
 * \brief Flow solvers and utilities.
 */

/**
 * \defgroup nssolver Navier-Stokes solver
 * \brief Implementation of parallel Navier-Stokes solver
 *
 * This is an example of using PetIBM to build a parallel incompressible flow
 * solver. The scheme used can be found in Perot (1993) and Chang et. al.
 * (2002). This example is a good starting point to learn how to use PetIBM to
 * build an immersed-boundary solver under Perot's framework.
 *
 * If readers are interested in using the Navier-Stokes solver instead of
 * coding, please refer to \ref md_doc_markdowns_runpetibm "Running PetIBM",
 * \ref md_doc_markdowns_examples2d "2D Examples", and
 * \ref md_doc_markdowns_examples3d "3D Examples".
 *
 * \b Reference: \n
 * \li Perot, J. B. (1993). An analysis of the fractional step method. Journal
 * of Computational Physics, 108(1), 51-58.
 * \li Chang, W., Giraldo, F., & Perot, B. (2002). Analysis of an exact
 * fractional step method. Journal of Computational Physics, 180(1), 183-199.
 *
 * \ingroup apps
 */

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    YAML::Node config;
    NavierStokesSolver solver;

    ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);
    ierr = PetscLogDefaultBegin(); CHKERRQ(ierr);

    // parse configuration files; store info in YAML node
    ierr = petibm::parser::getSettings(config); CHKERRQ(ierr);

    // initialize the Navier-Stokes solver
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

    // destroy the Navier-Stokes solver
    ierr = solver.destroy(); CHKERRQ(ierr);

    ierr = PetscFinalize(); CHKERRQ(ierr);

    return 0;
}  // main
