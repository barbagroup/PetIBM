/**
 * \file decoupledibpm/main.cpp
 * \brief Main function of decoupled IBPM solver (Li et. al. 2016).
 * \see decoupledibpm
 * \ingroup decoupledibpm
 */

#include <iomanip>
#include <sys/stat.h>

// PETSc
#include <petscsys.h>

// YAML-CPP
#include <yaml-cpp/yaml.h>

// PetIBM
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
 * \see nssolver, tairacolonius
 * \ingroup apps
 */


int main(int argc, char **argv)
{
    PetscErrorCode ierr;

    YAML::Node              config;
    petibm::type::Mesh      mesh;
    petibm::type::Boundary  bd;
    petibm::type::BodyPack  bodies;

    ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);
    ierr = PetscLogDefaultBegin(); CHKERRQ(ierr);
    
    // get all settings and save into `config`
    ierr = petibm::parser::getSettings(config); CHKERRQ(ierr);

    // create mesh
    ierr = petibm::mesh::createMesh(PETSC_COMM_WORLD, config, mesh); CHKERRQ(ierr);
    
    // output the data of the mesh
    ierr = mesh->write(config["directory"].as<std::string>() + "/grid"); CHKERRQ(ierr);
    
    // create bc
    ierr = petibm::boundary::createBoundary(mesh, config, bd); CHKERRQ(ierr);
    
    // create body pack
    ierr = petibm::body::createBodyPack(mesh, config, bodies); CHKERRQ(ierr);
    

    
    DecoupledIBPMSolver solver;

    // initialize the solver based on given mesh, boundary, bodies, and config
    ierr = solver.initialize(mesh, bd, bodies, config); CHKERRQ(ierr);

    

    // starting step
    PetscInt start = config["parameters"]["startStep"].as<PetscInt>();
    
    // end step
    PetscInt end = start + config["parameters"]["nt"].as<PetscInt>();
    
    // number of steps to save solutions
    PetscInt nsave = config["parameters"]["nsave"].as<PetscInt>();
    
    // number of steps to save solutions
    PetscInt nrestart = config["parameters"]["nrestart"].as<PetscInt>();

    // time-step size
    PetscReal dt = config["parameters"]["dt"].as<PetscReal>();

    // current time
    PetscReal t = 0.0;
             
    // directory where file to log information of linear solvers in
    std::string iterationsFile = config["directory"].as<std::string>() + "/";
             
    // directory where file to log averaged Lagrangian forces in
    std::string forceFile = config["directory"].as<std::string>() + "/";
    
    
    if (start == 0) // write initial solutions to a HDF5
    {
        ierr = PetscPrintf(PETSC_COMM_WORLD,
                "[time-step 0] Writing solution... "); CHKERRQ(ierr);
        
        std::stringstream ss;
        ss << "/" << std::setfill('0') << std::setw(7) << 0;
        ierr = solver.write(
                t, (config["solution"].as<std::string>() + ss.str()));
        CHKERRQ(ierr);
        
        ierr = PetscPrintf(PETSC_COMM_WORLD, "done\n"); CHKERRQ(ierr);

        iterationsFile += "iterations.txt";
        forceFile += "forces.txt";
    }
    else // restart
    {
        ierr = PetscPrintf(PETSC_COMM_WORLD,
                "[time-step %d] Read solution... ", start); CHKERRQ(ierr);
        
        std::stringstream ss;
        ss << "/" << std::setfill('0') << std::setw(7) << start;
        ierr = solver.readRestartData(
                (config["solution"].as<std::string>() + ss.str()), t);
        CHKERRQ(ierr);
        
        ierr = PetscPrintf(PETSC_COMM_WORLD, "done\n"); CHKERRQ(ierr);

        iterationsFile += "iterations-" + std::to_string(start) + ".txt";
        forceFile += "forces-" + std::to_string(start) + ".txt";
    }

    // initialize the PetscViewers of ASCII files in the solver class
    ierr = solver.initializeASCIIFiles(iterationsFile); CHKERRQ(ierr);
    ierr = solver.initializeASCIIFiles(forceFile); CHKERRQ(ierr);
    
    // start time marching
    for (int ite=start+1; ite<=end; ite++)
    {
        // update current time
        t += dt;
        
        // advance one time-step
        ierr = solver.advance(); CHKERRQ(ierr);
        ierr = solver.writeIterations(ite, iterationsFile); CHKERRQ(ierr);
        
        if (ite % nsave == 0)
        {
            ierr = PetscPrintf(PETSC_COMM_WORLD,
                    "[time-step %d] Writing solution... ", ite); CHKERRQ(ierr);
            
            std::stringstream ss;
            ss << "/" << std::setfill('0') << std::setw(7) << ite;
            ierr = solver.write(
                    t, (config["solution"].as<std::string>() + ss.str())); 
            CHKERRQ(ierr);
            
            ierr = PetscPrintf(PETSC_COMM_WORLD, "done\n"); CHKERRQ(ierr);
        }
        
        if (ite % nrestart == 0)
        {
            ierr = PetscPrintf(PETSC_COMM_WORLD, 
                    "[time-step %d] Writing necessary data for restarting... ",
                    ite); CHKERRQ(ierr);
            
            std::stringstream ss;
            ss << "/" << std::setfill('0') << std::setw(7) << ite;
            ierr = solver.writeRestartData(
                    t, (config["solution"].as<std::string>() + ss.str()));
            CHKERRQ(ierr);
            
            ierr = PetscPrintf(PETSC_COMM_WORLD, "done\n"); CHKERRQ(ierr);
        }
        
        // write averaged force
        ierr = solver.writeIntegratedForces(t, forceFile); CHKERRQ(ierr);
    }

    ierr = PetscFinalize(); CHKERRQ(ierr);
    
    return 0;
} // main
