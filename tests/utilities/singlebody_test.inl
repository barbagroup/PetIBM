/**
 * \file singlebody_test.cpp
 * \brief Unit-tests for the class `SingleBody`.
 */

#include <vector>

#include "gtest/gtest.h"
#include "yaml-cpp/yaml.h"

#include "utilities/SingleBody.h"
#include "utilities/parser.h"
#include "utilities/FlowDescription.h"

using namespace petibm::utilities;


class SingleBodyTest : public ::testing::Test
{
protected:
	
	SingleBodyTest(){  };
	
	virtual ~SingleBodyTest(){  };

	virtual void SetUp()
	{
		YAML::Node config;
		parser::parseYAMLConfigFile("data/config.yaml", config);
		FlowDescription flow = FlowDescription(PETSC_COMM_WORLD,
		                                       config["flowDescription"]);
		mesh = CartesianMesh(PETSC_COMM_WORLD,
		                     config["cartesianMesh"],
		                     flow.BCInfo,
		                     types::str2out["Binary"]);
		body = SingleBody(mesh, "data/body.txt", "body");
	};
	
	virtual void TearDown(){  };

	SingleBody body;
	CartesianMesh mesh;

}; // SingleBodyTest


TEST_F(SingleBodyTest, initWithFilePath2D)
{
	body.init(mesh, "data/body.txt", "body");
	ASSERT_EQ(2, body.dim);
	ASSERT_EQ("body", body.name);
	ASSERT_EQ(4, body.nPts);
	std::vector<PetscReal> xCoords = {0.25, 0.75, 0.75, 0.25},
	                       yCoords = {0.25, 0.25, 0.75, 0.75};
	for (unsigned int i=0; i<4; i++)
	{
		ASSERT_EQ(xCoords[i], body.coords[i][0]);
		ASSERT_EQ(yCoords[i], body.coords[i][1]);
	}
	PetscMPIInt size;
	MPI_Comm_size(PETSC_COMM_WORLD, &size);
	if (size == 1)
		ASSERT_EQ(body.nPts, body.nLclPts);
}


TEST_F(SingleBodyTest, findProc)
{
	PetscMPIInt size, index;
	MPI_Comm_size(PETSC_COMM_WORLD, &size);
	if (size == 1)
		for (int i=0; i<body.nPts; i++)
		{
			body.findProc(i, index);
			ASSERT_EQ(0, index);
		}
}


TEST_F(SingleBodyTest, getGlobalIndex)
{
	PetscInt globalIndex, counter = 0;
	PetscInt ndof = 2;
	for (int i=0; i<body.nPts; i++)
		for (int d=0; d<ndof; d++)
		{
			body.getGlobalIndex(i, d, globalIndex);
			ASSERT_EQ(counter, globalIndex);
			MatStencil stencil = {0, 0, i, d};
			body.getGlobalIndex(stencil, globalIndex);
			ASSERT_EQ(counter, globalIndex);
			counter++;
		}
}


TEST_F(SingleBodyTest, calculateAvgForces)
{
	Vec f;
	DMCreateGlobalVector(body.da, &f);
	VecSet(f, 1.0);
	types::RealVec1D avg(body.dim);
	body.calculateAvgForces(f, avg);
	ASSERT_EQ(body.dim, (int)avg.size());
	for (unsigned int d=0; d<avg.size(); d++)
		ASSERT_EQ(body.nPts * 1.0, avg[d]);
	VecDestroy(&f);
}
