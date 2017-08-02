/**
 * \file singlebody_test.cpp
 * \brief Unit-tests for the class `SingleBody`.
 */

#include <vector>

#include <petsc.h>

#include "gtest/gtest.h"
#include "yaml-cpp/yaml.h"

#include "utilities/SingleBody.h"
#include "utilities/parser.h"
#include "utilities/FlowDescription.h"
#include "utilities/CartesianMesh.h"

using namespace petibm::utilities;


class SingleBodyTest : public ::testing::Test
{
protected:
	
	SingleBodyTest(){  };
	
	virtual ~SingleBodyTest(){  };

	virtual void SetUp()
	{
		YAML::Node config;
		FlowDescription flow;
		// create 2D mesh and body
		parser::parseYAMLConfigFile("data/config2d.yaml", config);
		flow = FlowDescription(PETSC_COMM_WORLD, config["flowDescription"]);
		mesh2d = CartesianMesh(PETSC_COMM_WORLD,
		                       config["cartesianMesh"],
		                       flow.BCInfo,
		                       types::str2out["Binary"]);
		body2d = SingleBody(mesh2d, "data/body2d.txt", "body2d");
		// create 3D mesh and body
		parser::parseYAMLConfigFile("data/config3d.yaml", config);
		flow = FlowDescription(PETSC_COMM_WORLD, config["flowDescription"]);
		mesh3d = CartesianMesh(PETSC_COMM_WORLD,
		                       config["cartesianMesh"],
		                       flow.BCInfo,
		                       types::str2out["Binary"]);
		body3d = SingleBody(mesh3d, "data/body3d.txt", "body3d");
	};
	
	virtual void TearDown(){  };

	SingleBody body2d, body3d;
	CartesianMesh mesh2d, mesh3d;

}; // SingleBodyTest


TEST_F(SingleBodyTest, initWithFilePath2D)
{
	body2d.init(mesh2d, "data/body2d.txt", "body2d");
	ASSERT_EQ(2, body2d.dim);
	ASSERT_EQ("body2d", body2d.name);
	ASSERT_EQ(4, body2d.nPts);
	std::vector<PetscReal> xCoords = {0.25, 0.75, 0.75, 0.25},
	                       yCoords = {0.25, 0.25, 0.75, 0.75};
	for (unsigned int i=0; i<4; i++)
	{
		ASSERT_EQ(xCoords[i], body2d.coords[i][0]);
		ASSERT_EQ(yCoords[i], body2d.coords[i][1]);
	}
	PetscMPIInt size;
	MPI_Comm_size(PETSC_COMM_WORLD, &size);
	if (size == 1)
		ASSERT_EQ(body2d.nPts, body2d.nLclPts);
}


TEST_F(SingleBodyTest, initWithFilePath3D)
{
	body3d.init(mesh3d, "data/body3d.txt", "body3d");
	ASSERT_EQ(3, body3d.dim);
	ASSERT_EQ("body3d", body3d.name);
	ASSERT_EQ(8, body3d.nPts);
	std::vector<PetscReal> xCoords = {0.25, 0.75, 0.75, 0.25, 0.25, 0.75, 0.75, 0.25},
	                       yCoords = {0.25, 0.25, 0.75, 0.75, 0.25, 0.25, 0.75, 0.75},
	                       zCoords = {0.25, 0.25, 0.25, 0.25, 0.75, 0.75, 0.75, 0.75};
	for (unsigned int i=0; i<4; i++)
	{
		ASSERT_EQ(xCoords[i], body3d.coords[i][0]);
		ASSERT_EQ(yCoords[i], body3d.coords[i][1]);
		ASSERT_EQ(zCoords[i], body3d.coords[i][2]);
	}
	PetscMPIInt size;
	MPI_Comm_size(PETSC_COMM_WORLD, &size);
	if (size == 1)
		ASSERT_EQ(body3d.nPts, body3d.nLclPts);
}


TEST_F(SingleBodyTest, findProc2D)
{
	PetscMPIInt size, index;
	MPI_Comm_size(PETSC_COMM_WORLD, &size);
	if (size == 1)
		for (int i=0; i<body2d.nPts; i++)
		{
			body2d.findProc(i, index);
			ASSERT_EQ(0, index);
		}
}


TEST_F(SingleBodyTest, findProc3D)
{
	PetscMPIInt size, index;
	MPI_Comm_size(PETSC_COMM_WORLD, &size);
	if (size == 1)
		for (int i=0; i<body3d.nPts; i++)
		{
			body3d.findProc(i, index);
			ASSERT_EQ(0, index);
		}
}


TEST_F(SingleBodyTest, getGlobalIndex2D)
{
	PetscInt globalIndex, counter = 0;
	PetscInt ndof = 2;
	for (int i=0; i<body2d.nPts; i++)
		for (int d=0; d<ndof; d++)
		{
			body2d.getGlobalIndex(i, d, globalIndex);
			ASSERT_EQ(counter, globalIndex);
			MatStencil stencil = {0, 0, i, d};
			body2d.getGlobalIndex(stencil, globalIndex);
			ASSERT_EQ(counter, globalIndex);
			counter++;
		}
}


TEST_F(SingleBodyTest, getGlobalIndex3D)
{
	PetscInt globalIndex, counter = 0;
	PetscInt ndof = 3;
	for (int i=0; i<body3d.nPts; i++)
		for (int d=0; d<ndof; d++)
		{
			body3d.getGlobalIndex(i, d, globalIndex);
			ASSERT_EQ(counter, globalIndex);
			MatStencil stencil = {0, 0, i, d};
			body3d.getGlobalIndex(stencil, globalIndex);
			ASSERT_EQ(counter, globalIndex);
			counter++;
		}
}


TEST_F(SingleBodyTest, calculateAvgForces2D)
{
	Vec f;
	DMCreateGlobalVector(body2d.da, &f);
	VecSet(f, 1.0);
	types::RealVec1D avg(2);
	body2d.calculateAvgForces(f, avg);
	ASSERT_EQ(body2d.dim, (int)avg.size());
	for (unsigned int d=0; d<avg.size(); d++)
		ASSERT_EQ(body2d.nPts * 1.0, avg[d]);
	VecDestroy(&f);
}


TEST_F(SingleBodyTest, calculateAvgForces3D)
{
	Vec f;
	DMCreateGlobalVector(body3d.da, &f);
	VecSet(f, 1.0);
	types::RealVec1D avg(3);
	body3d.calculateAvgForces(f, avg);
	ASSERT_EQ(body3d.dim, (int)avg.size());
	for (unsigned int d=0; d<avg.size(); d++)
		ASSERT_EQ(body3d.nPts * 1.0, avg[d]);
	VecDestroy(&f);
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
