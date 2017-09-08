/**
 * \file cartesianmesh_test.cpp
 * \brief Unit-tests for the class `CartesianMesh`.
 */

#include <petsc.h>

#include "gtest/gtest.h"
#include "yaml-cpp/yaml.h"

#include "petibm/cartesianmesh.h"
#include "petibm/parser.h"
#include "petibm/flowdescription.h"

using namespace petibm::utilities;


class CartesianMeshTest : public ::testing::Test
{
protected:

	CartesianMeshTest(){  };

	virtual ~CartesianMeshTest(){  };

	virtual void SetUp()
	{
		YAML::Node config;
		FlowDescription flow;
		parser::parseYAMLConfigFile("utilities/data/config3d.yaml", config);
		flow = FlowDescription(PETSC_COMM_WORLD, config["flowDescription"]);
		mesh = CartesianMesh(PETSC_COMM_WORLD,
		                     config["cartesianMesh"],
		                     flow.BCInfo);

	};

	virtual void TearDown(){  };

	CartesianMesh mesh;
	std::string directory = ".";
}; // CartesianMeshTest


TEST_F(CartesianMeshTest, init3D)
{
	YAML::Node config;
	FlowDescription flow;
	parser::parseYAMLConfigFile("utilities/data/config3d.yaml", config);
	flow = FlowDescription(PETSC_COMM_WORLD, config["flowDescription"]);
	mesh.init(PETSC_COMM_WORLD,
            config["cartesianMesh"],
            flow.BCInfo);
	// check dimension
	ASSERT_EQ(3, mesh.dim);
	// check domain limits
	for (unsigned int d=0; d<3; d++)
	{
		ASSERT_EQ(0.1, mesh.min[d]);
		ASSERT_EQ(1.1, mesh.max[d]);
	}
	// check number of points for each field
	std::vector<PetscReal> n_exp(3);
	// for x-velocity
	n_exp = {9, 10, 10};
	for (unsigned int d=0; d<3; d++)
		ASSERT_EQ(n_exp[d], mesh.n[0][d]);
	// for y-velocity
	n_exp = {10, 9, 10};
	for (unsigned int d=0; d<3; d++)
		ASSERT_EQ(n_exp[d], mesh.n[1][d]);
	// for z-velocity
	n_exp = {10, 10, 10};
	for (unsigned int d=0; d<3; d++)
		ASSERT_EQ(n_exp[d], mesh.n[2][d]);
	// for pressure
	n_exp = {10, 10, 10};
	for (unsigned int d=0; d<3; d++)
		ASSERT_EQ(n_exp[d], mesh.n[3][d]);
	// for vertices
	n_exp = {11, 11, 11};
	for (unsigned int d=0; d<3; d++)
		ASSERT_EQ(n_exp[d], mesh.n[4][d]);
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
