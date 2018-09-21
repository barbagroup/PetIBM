/**
 * \file decoupledibpm/main.cpp
 * \brief Main function of decoupled IBPM solver (Li et. al. 2016).
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 * \see decoupledibpm
 * \ingroup decoupledibpm
 */

#include <sys/stat.h>
#include <iomanip>

#include <petscsys.h>
#include <yaml-cpp/yaml.h>

#include <petibm/io.h>
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
    petibm::type::Mesh mesh;
    petibm::type::Boundary bd;
    petibm::type::BodyPack bodies;
    DecoupledIBPMSolver solver;
    std::string filePath, iterationsFile, forcesFile;

    ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);
    ierr = PetscLogDefaultBegin(); CHKERRQ(ierr);

    // parse configuration files; store info in YAML node
    ierr = petibm::parser::getSettings(config); CHKERRQ(ierr);

    // create mesh
    ierr = petibm::mesh::createMesh(PETSC_COMM_WORLD, config, mesh);
    CHKERRQ(ierr);

    // write grid data into HDF5 file
    filePath = config["directory"].as<std::string>() + "/grid.h5";
    ierr = mesh->write(filePath); CHKERRQ(ierr);

    // create bc
    ierr = petibm::boundary::createBoundary(mesh, config, bd); CHKERRQ(ierr);

    // create body pack
    ierr = petibm::body::createBodyPack(PETSC_COMM_WORLD, mesh->dim, config,
                                        bodies); CHKERRQ(ierr);
    ierr = bodies->updateMeshIdx(mesh); CHKERRQ(ierr);

    // initialize the solver
    ierr = solver.initialize(mesh, bd, bodies, config); CHKERRQ(ierr);

    // store the value of some temporal parameters for readability
    // starting time-step index
    PetscInt nstart = config["parameters"]["startStep"].as<PetscInt>();
    // number of time steps to compute
    PetscInt nt = config["parameters"]["nt"].as<PetscInt>();
    // frequency at which the flow solution is saved (in number of time steps)
    PetscInt nsave = config["parameters"]["nsave"].as<PetscInt>();
    // frequency at which necessary files to restart a run are saved
    //(in number of time steps)
    PetscInt nrestart = config["parameters"]["nrestart"].as<PetscInt>();
    // time-step size
    PetscReal dt = config["parameters"]["dt"].as<PetscReal>();
    // current time set to initial time
    PetscReal t = config["parameters"]["t"].as<PetscReal>(0.0);

    if (nstart == 0)  // write initial solution into HDF5 file
    {
        ierr =
            PetscPrintf(PETSC_COMM_WORLD, "[time-step 0] Writing solution... ");
        CHKERRQ(ierr);

        std::stringstream ss;
        ss << std::setfill('0') << std::setw(7) << 0;
        filePath =
            config["solution"].as<std::string>() + "/" + ss.str() + ".h5";
        ierr = solver.write(t, filePath); CHKERRQ(ierr);

        ierr = PetscPrintf(PETSC_COMM_WORLD, "done\n"); CHKERRQ(ierr);

        iterationsFile =
            config["directory"].as<std::string>() + "/iterations.txt";
        forcesFile = config["directory"].as<std::string>() + "/forces.txt";
    }
    else  // read solution and restart variables from HDF5 file
    {
        ierr = PetscPrintf(PETSC_COMM_WORLD, "[time-step %d] Read solution... ",
                           nstart); CHKERRQ(ierr);

        std::stringstream ss;
        ss << std::setfill('0') << std::setw(7) << nstart;
        filePath =
            config["solution"].as<std::string>() + "/" + ss.str() + ".h5";
        ierr = solver.readRestartData(filePath, t); CHKERRQ(ierr);

        ierr = PetscPrintf(PETSC_COMM_WORLD, "done\n"); CHKERRQ(ierr);

        iterationsFile = config["directory"].as<std::string>() +
                         "/iterations-" + std::to_string(nstart) + ".txt";
        forcesFile = config["directory"].as<std::string>() + "/forces-" +
                     std::to_string(nstart) + ".txt";
    }

    // initialize the PetscViewer objects for ASCII files
    ierr = solver.initializeASCIIFiles(iterationsFile); CHKERRQ(ierr);
    ierr = solver.initializeASCIIFiles(forcesFile); CHKERRQ(ierr);

    // Time integration
    for (int ite = nstart + 1; ite <= nstart + nt; ite++)
    {
        // advance solution
        t += dt;
        ierr = solver.advance(); CHKERRQ(ierr);
        // write solvers iterations into ASCII file
        ierr = solver.writeIterations(ite, iterationsFile); CHKERRQ(ierr);
        // write hydrodynamic forces into ASCII file
        ierr = solver.writeIntegratedForces(t, forcesFile); CHKERRQ(ierr);

        if (ite % nsave == 0)  // write flow solution into HDF5 file
        {
            std::stringstream ss;
            ss << std::setfill('0') << std::setw(7) << ite;
            filePath =
                config["solution"].as<std::string>() + "/" + ss.str() + ".h5";

            ierr = PetscPrintf(PETSC_COMM_WORLD,
                               "[time-step %d] Writing solution... ", ite);
            CHKERRQ(ierr);
            ierr = solver.write(t, filePath); CHKERRQ(ierr);
            ierr = PetscPrintf(PETSC_COMM_WORLD, "done\n"); CHKERRQ(ierr);

            if (ite % nrestart == 0)  // write restart data into HDF5 file
            {
                ierr =
                    PetscPrintf(PETSC_COMM_WORLD,
                                "[time-step %d] Writing restart data... ", ite);
                CHKERRQ(ierr);
                ierr = solver.writeRestartData(t, filePath); CHKERRQ(ierr);
                ierr = PetscPrintf(PETSC_COMM_WORLD, "done\n"); CHKERRQ(ierr);
            }

            // write summary of PETSc logging into ASCII file
            filePath =
                config["solution"].as<std::string>() + "/" + ss.str() + ".log";
            ierr = petibm::io::writePetscLog(PETSC_COMM_WORLD, filePath);
            CHKERRQ(ierr);
        }
    }

    // manually destroy PETSc objects before PetscFinalize
    ierr = solver.destroy(); CHKERRQ(ierr);
    ierr = bodies->destroy(); CHKERRQ(ierr);
    ierr = bd->destroy(); CHKERRQ(ierr);
    ierr = mesh->destroy(); CHKERRQ(ierr);

    ierr = PetscFinalize(); CHKERRQ(ierr);

    return 0;
}  // main
