/**
 * \file navierstokes/main.cpp
 * \brief Main function of the Navier-Stokes solver.
 * \see nssolver
 * \ingroup nssolver
 */

#include <iomanip>
#include <sys/stat.h>

// PETSc
#include <petscsys.h>

// YAML-CPP
#include <yaml-cpp/yaml.h>

// PetIBM
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
 * solver. The scheme used can be found in Perot (1993) and Chang et. al. (2002).
 * This example is a good starting point to learn how to use PetIBM to build
 * an immersed-boundary solver under Perot's framework.
 * 
 * If readers are interested in using the Navier-Stokes solver instead of coding,
 * please refer to 
 * \ref md_doc_markdowns_runpetibm "Running PetIBM",
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

    YAML::Node              config;
    petibm::type::Mesh      mesh;
    petibm::type::Boundary  bd;

    ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);
    
    // get all settings and save into `config`
    ierr = petibm::parser::getSettings(config); CHKERRQ(ierr);

    // create mesh
    ierr = petibm::mesh::createMesh(PETSC_COMM_WORLD, config, mesh); CHKERRQ(ierr);
    
    // output the data of the mesh
    ierr = mesh->write(config["directory"].as<std::string>() + "/grid"); CHKERRQ(ierr);
    
    // create bc
    ierr = petibm::boundary::createBoundary(mesh, config, bd); CHKERRQ(ierr);

    
    
    NavierStokesSolver      solver;
    
    // initialize the flow solver based on given mesh, boundary, and config
    ierr = solver.initialize(mesh, bd, config); CHKERRQ(ierr);

    

    // starting step
    PetscInt start = config["parameters"]["startStep"].as<PetscInt>();
    
    // end step
    PetscInt end = start + config["parameters"]["nt"].as<PetscInt>();
    
    // number of steps to save solutions
    PetscInt nsave = config["parameters"]["nsave"].as<PetscInt>();
    
    // number of steps to save solutions
    PetscInt nrestart = config["parameters"]["nrestart"].as<PetscInt>();
             
    // directory where file to log information of linear solvers in
    std::string iterationsFile = config["directory"].as<std::string>() + "/";
    
    
    if (start == 0) // write initial solutions to a HDF5
    {
        ierr = PetscPrintf(PETSC_COMM_WORLD,
                "[time-step 0] Writing solution... "); CHKERRQ(ierr);
        
        std::stringstream ss;
        ss << "/" << std::setfill('0') << std::setw(7) << 0;
        ierr = solver.write(
                (config["solution"].as<std::string>() + ss.str()));
        CHKERRQ(ierr);
        
        ierr = PetscPrintf(PETSC_COMM_WORLD, "done\n"); CHKERRQ(ierr);

             
        iterationsFile += "iterations.txt";
    }
    else // restart
    {
        ierr = PetscPrintf(PETSC_COMM_WORLD,
                "[time-step %d] Read solution... ", start); CHKERRQ(ierr);
        
        std::stringstream ss;
        ss << "/" << std::setfill('0') << std::setw(7) << start;
        ierr = solver.readRestartData(
                (config["solution"].as<std::string>() + ss.str()));
        CHKERRQ(ierr);
        
        ierr = PetscPrintf(PETSC_COMM_WORLD, "done\n"); CHKERRQ(ierr);

        iterationsFile += "iterations-" + std::to_string(start) + ".txt";
    }

    // initialize the PetscViewers of ASCII files in the solver class
    ierr = solver.initializeASCIIFiles(iterationsFile); CHKERRQ(ierr);
    
    // start time marching
    for (int ite=start+1; ite<=end; ite++)
    {
        ierr = solver.advance(); CHKERRQ(ierr);
        ierr = solver.writeIterations(ite, iterationsFile); CHKERRQ(ierr);
        
        if (ite % nsave == 0)
        {
            ierr = PetscPrintf(PETSC_COMM_WORLD, 
                    "[time-step %d] Writing solution... ", ite); CHKERRQ(ierr);
            
            std::stringstream ss;
            ss << "/" << std::setfill('0') << std::setw(7) << ite;
            ierr = solver.write(
                    (config["solution"].as<std::string>() + ss.str()));
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
                    (config["solution"].as<std::string>() + ss.str()));
            CHKERRQ(ierr);
            
            ierr = PetscPrintf(PETSC_COMM_WORLD, "done\n"); CHKERRQ(ierr);
        }
    }

    ierr = PetscFinalize(); CHKERRQ(ierr);
    
    return 0;
} // main
