/**
 * \file singleboundary_test.cpp
 * \brief Unit-tests for the class `SingleBoundary`.
 */

#include <petsc.h>

#include "gtest/gtest.h"
#include "yaml-cpp/yaml.h"

#include "utilities/SingleBoundary.h"
#include "utilities/parser.h"
#include "utilities/FlowDescription.h"
#include "utilities/CartesianMesh.h"

using namespace petibm::utilities;


class SingleBoundaryTest : public ::testing::Test
{
protected:

	SingleBoundaryTest(){  };

	virtual ~SingleBoundaryTest(){  };

	virtual void SetUp()
	{
		YAML::Node config;
		FlowDescription flow;
		CartesianMesh mesh;
		parser::parseYAMLConfigFile("data/config3d.yaml", config);
		flow = FlowDescription(PETSC_COMM_WORLD, config["flowDescription"]);
		mesh = CartesianMesh(PETSC_COMM_WORLD,
		                     config["cartesianMesh"],
		                     flow.BCInfo);
		boundary = SingleBoundary(mesh, types::str2bl["XMINUS"]);
	};

	virtual void TearDown(){  };

	SingleBoundary boundary;
}; // SingleBoundaryTest


TEST_F(SingleBoundaryTest, init)
{
	YAML::Node config;
	parser::parseYAMLConfigFile("data/config3d.yaml", config);
	FlowDescription flow = FlowDescription(PETSC_COMM_WORLD,
	                                       config["flowDescription"]);
	CartesianMesh mesh = CartesianMesh(PETSC_COMM_WORLD,
								                     config["cartesianMesh"],
								                     flow.BCInfo);
	boundary.init(mesh, types::BCLoc::XMINUS);
	ASSERT_EQ(3, boundary.dim);
	ASSERT_EQ(types::BCLoc::XMINUS, boundary.loc);
	PetscMPIInt size;
	MPI_Comm_size(PETSC_COMM_WORLD, &size);
	if (size == 1)
		ASSERT_EQ(PETSC_TRUE, boundary.onThisProc);
	ASSERT_EQ(3, (int)boundary.type.size());
	for (unsigned int d=0; d<3; d++)
	{
		ASSERT_EQ(types::BCType::DIRICHLET, boundary.type[d]);
		ASSERT_EQ(0.0, boundary.value[d]);
	}
}


// Run all tests
int main(int argc, char **argv)
{
	PetscErrorCode ierr, status;

	::testing::InitGoogleTest(&argc, argv);
	ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);
	status = RUN_ALL_TESTS();
	ierr = PetscFinalize(); CHKERRQ(ierr);

	return status;
} // main
