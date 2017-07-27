/** Decoupled IBPM solver.
 * \file main.cpp
 */

#include <petsc.h>

#include <yaml-cpp/yaml.h>

#include "decoupledibpm.h"
#include "CartesianMesh.h"
#include "FlowDescription.h"
#include "SimulationParameters.h"
#include "parser.h"


int main(int argc, char **argv)
{
	PetscErrorCode ierr;

	ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);

	SimulationParameters params;
	FlowDescription flow;
	CartesianMesh mesh;
	YAML::Node config;
	char path[PETSC_MAX_PATH_LEN];
	std::string directory,
							configpath;
	PetscBool flg;

	// parse command-line looking for simulation directory and configuration path
	directory = ".";
	ierr = PetscOptionsGetString(nullptr, nullptr, "-directory",
	                             path, sizeof(path), &flg); CHKERRQ(ierr);
	if (flg)
		directory = path;
	configpath = directory + "/config.yaml";
	ierr = PetscOptionsGetString(nullptr, nullptr,"-config",
	                             path, sizeof(path), &flg); CHKERRQ(ierr);
	if (flg)
		configpath = path;

	// parse configuration file
	ierr = parser::parseYAMLConfigFile(configpath, config); CHKERRQ(ierr);

	ierr = params.init(PETSC_COMM_WORLD,
	                   config["simulationParameters"], directory); CHKERRQ(ierr);
	ierr = params.printInfo(); CHKERRQ(ierr);

	ierr = flow.init(PETSC_COMM_WORLD,
	                 config["flowDescription"]); CHKERRQ(ierr);
	ierr = flow.printInfo(); CHKERRQ(ierr);

	ierr = mesh.init(PETSC_COMM_WORLD,
	                 config["cartesianMesh"], flow.BCInfo,
	                 params.output.format); CHKERRQ(ierr);
	ierr = mesh.printInfo(); CHKERRQ(ierr);
	ierr = mesh.write(params.caseDir, "grid"); CHKERRQ(ierr);
	
	DecoupledIBPMSolver solver = DecoupledIBPMSolver(
			mesh, flow, params); CHKERRQ(ierr);

	ierr = solver.initialize(); CHKERRQ(ierr);

	PetscInt start = params.step.nStart,
	         end = params.step.nStart + params.step.nTotal,
	         nsave = params.step.nSave;
	std::string iterationsFile = directory + "/iterations.txt";
	for (int ite=start+1; ite<=end; ite++)
	{
		ierr = solver.solve(); CHKERRQ(ierr);
		if (ite % nsave == 0)
		{
			ierr = solver.write(
					config["caseDir"].as<std::string>(), "01"); CHKERRQ(ierr);
		}
		ierr = solver.writeIterations(ite, iterationsFile); CHKERRQ(ierr);
	}

	ierr = solver.finalize(); CHKERRQ(ierr);

	ierr = PetscFinalize(); CHKERRQ(ierr);
	
	return 0;
} // main
