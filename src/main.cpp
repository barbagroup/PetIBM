#ifndef DIMENSIONS
#define DIMENSIONS 2
#endif

#include <petsc.h>
#include <yaml-cpp/yaml.h>

#include "navierstokes.h"
#include "CartesianMesh.h"
#include "FlowDescription.h"
#include "SimulationParameters.h"
#include "parser.h"


int main(int argc, char **argv)
{
	PetscErrorCode ierr;

	ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);

	YAML::Node config;
	ierr = parser::parseYAMLConfigFile("./yamls/lidDriven/2d",
	                                   config); CHKERRQ(ierr);

	SimulationParameters params;
	ierr = params.init(PETSC_COMM_WORLD, config["simulationParameters"],
	                   config["caseDir"].as<std::string>()); CHKERRQ(ierr);
	ierr = params.printInfo(); CHKERRQ(ierr);

	FlowDescription flow;
	ierr = flow.init(PETSC_COMM_WORLD, config["flowDescription"]); CHKERRQ(ierr);
	ierr = flow.printInfo(); CHKERRQ(ierr);

	CartesianMesh mesh;
	ierr = mesh.init(PETSC_COMM_WORLD, config["cartesianMesh"], flow.BCInfo,
	                 params.output.format); CHKERRQ(ierr);
	ierr = mesh.printInfo(); CHKERRQ(ierr);
	ierr = mesh.write(params.caseDir, "grid"); CHKERRQ(ierr);
	
	NavierStokesSolver solver = 
			NavierStokesSolver(mesh, flow, params); CHKERRQ(ierr);

	ierr = solver.initialize(); CHKERRQ(ierr);

	PetscInt start = params.step.nStart,
					 end = params.step.nStart + params.step.nTotal,
					 nsave = params.step.nSave;
	std::string iterationsFile = 
			config["caseDir"].as<std::string>() + "/iterations.txt";
	for (int ite=start; ite<end; ite++)
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
