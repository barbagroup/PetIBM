/** Navier-Stokes solver
 * \file main.cpp
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
    // TODO: ierr = mesh.write(params.caseDir, "grid"); CHKERRQ(ierr);
    
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
             
    // file to log information of linear solvers
    std::string iterationsFile = 
        config["directory"].as<std::string>() + "/iterations.txt";
    
    
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
    }
    
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
