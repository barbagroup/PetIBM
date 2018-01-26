/**
 * \file singlebody_test.cpp
 * \brief Unit-tests for the class `SingleBody`.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 */

#include <vector>

#include <petsc.h>

#include <gtest/gtest.h>
#include <yaml-cpp/yaml.h>

#include <petibm/mesh.h>
#include <petibm/singlebody.h>

using namespace petibm;


class SingleBodyTest : public ::testing::Test
{
protected:
    
    SingleBodyTest(){  };
    
    virtual ~SingleBodyTest(){  };

    virtual void SetUp()
    {
        using namespace YAML;
        YAML::Node config;
        
        config["mesh"].push_back(Node(NodeType::Map));
        config["mesh"][0]["direction"] = "x";
        config["mesh"][1]["direction"] = "y";
        for(unsigned int i=0; i<2; ++i)
        {
            config["mesh"][i]["start"] = "0.1";
            config["mesh"][i]["subDomains"].push_back(Node(NodeType::Map));
            config["mesh"][i]["subDomains"][0]["end"] = 1.0;
            config["mesh"][i]["subDomains"][0]["cells"] = 10;
            config["mesh"][i]["subDomains"][0]["stretchRatio"] = 1.0;
        }
        
        config["flow"] = YAML::Node(NodeType::Map);
        config["flow"]["boundaryConditions"].push_back(Node(NodeType::Map));
        config["flow"]["boundaryConditions"][0]["location"] = "xMinus";
        config["flow"]["boundaryConditions"][1]["location"] = "xPlus";
        config["flow"]["boundaryConditions"][2]["location"] = "yMinus";
        config["flow"]["boundaryConditions"][3]["location"] = "yPlus";
        
        for(unsigned int i=0; i<4; ++i)
        {
            config["flow"]["boundaryConditions"][i]["u"][0] = "DIRICHLET";
            config["flow"]["boundaryConditions"][i]["u"][1] = 0.0;
            config["flow"]["boundaryConditions"][i]["v"][0] = "DIRICHLET";
            config["flow"]["boundaryConditions"][i]["v"][1] = 0.0;
        }
        config["flow"]["boundaryConditions"][3]["u"][1] = 1.0;
        
        config["directory"] = "body/";
        config["bodies"][0]["type"] = "points";
        config["bodies"][0]["file"] = "body2d.txt";
        
        
        // create 2D mesh and body
        petibm::mesh::createMesh(PETSC_COMM_WORLD, config, mesh2d);
        petibm::body::createSingleBody(mesh2d, "points", "body2d01", "body/body2d.txt", body2d);
        
        // create 3D mesh and body
        config["mesh"][2]["direction"] = "z";
        config["mesh"][2]["start"] = "0.1";
        config["mesh"][2]["subDomains"].push_back(Node(NodeType::Map));
        config["mesh"][2]["subDomains"][0]["end"] = 1.0;
        config["mesh"][2]["subDomains"][0]["cells"] = 10;
        config["mesh"][2]["subDomains"][0]["stretchRatio"] = 1.0;
        for(unsigned int i=0; i<4; ++i)
        {
            config["flow"]["boundaryConditions"][i]["w"][0] = "DIRICHLET";
            config["flow"]["boundaryConditions"][i]["w"][1] = 0.0;
        }
        config["flow"]["boundaryConditions"][4]["location"] = "zMinus";
        config["flow"]["boundaryConditions"][5]["location"] = "zPlus";
        
        for(unsigned int i=4; i<6; ++i)
        {
            config["flow"]["boundaryConditions"][i]["u"][0] = "PERIODIC";
            config["flow"]["boundaryConditions"][i]["u"][1] = 0.0;
            config["flow"]["boundaryConditions"][i]["v"][0] = "PERIODIC";
            config["flow"]["boundaryConditions"][i]["v"][1] = 0.0;
            config["flow"]["boundaryConditions"][i]["w"][0] = "PERIODIC";
            config["flow"]["boundaryConditions"][i]["w"][1] = 0.0;
        }
        config["bodies"][0]["file"] = "body3d.txt";
        
        petibm::mesh::createMesh(PETSC_COMM_WORLD, config, mesh3d);
        petibm::body::createSingleBody(mesh3d, "points", "body3d01", "body/body3d.txt", body3d);
    };
    
    virtual void TearDown(){  };

    type::SingleBody body2d, body3d;
    type::Mesh mesh2d, mesh3d;

}; // SingleBodyTest


TEST_F(SingleBodyTest, initWithFilePath2D)
{
    ASSERT_EQ(2, body2d->dim);
    ASSERT_EQ("body2d01", body2d->name);
    ASSERT_EQ(4, body2d->nPts);
    std::vector<PetscReal> xCoords = {0.25, 0.75, 0.75, 0.25},
                           yCoords = {0.25, 0.25, 0.75, 0.75};
    for (unsigned int i=0; i<4; i++)
    {
        ASSERT_EQ(xCoords[i], body2d->coords[i][0]);
        ASSERT_EQ(yCoords[i], body2d->coords[i][1]);
    }
    PetscMPIInt size;
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    if (size == 1)
        ASSERT_EQ(body2d->nPts, body2d->nLclPts);
}


TEST_F(SingleBodyTest, initWithFilePath3D)
{
    ASSERT_EQ(3, body3d->dim);
    ASSERT_EQ("body3d01", body3d->name);
    ASSERT_EQ(8, body3d->nPts);
    std::vector<PetscReal> xCoords = {0.25, 0.75, 0.75, 0.25, 0.25, 0.75, 0.75, 0.25},
                           yCoords = {0.25, 0.25, 0.75, 0.75, 0.25, 0.25, 0.75, 0.75},
                           zCoords = {0.25, 0.25, 0.25, 0.25, 0.75, 0.75, 0.75, 0.75};
    for (unsigned int i=0; i<4; i++)
    {
        ASSERT_EQ(xCoords[i], body3d->coords[i][0]);
        ASSERT_EQ(yCoords[i], body3d->coords[i][1]);
        ASSERT_EQ(zCoords[i], body3d->coords[i][2]);
    }
    PetscMPIInt size;
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    if (size == 1)
        ASSERT_EQ(body3d->nPts, body3d->nLclPts);
}


TEST_F(SingleBodyTest, findProc2D)
{
    PetscMPIInt size, index;
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    if (size == 1)
        for (int i=0; i<body2d->nPts; i++)
        {
            body2d->findProc(i, index);
            ASSERT_EQ(0, index);
        }
}


TEST_F(SingleBodyTest, findProc3D)
{
    PetscMPIInt size, index;
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    if (size == 1)
        for (int i=0; i<body3d->nPts; i++)
        {
            body3d->findProc(i, index);
            ASSERT_EQ(0, index);
        }
}


TEST_F(SingleBodyTest, getGlobalIndex2D)
{
    PetscInt globalIndex, counter = 0;
    PetscInt ndof = 2;
    for (int i=0; i<body2d->nPts; i++)
        for (int d=0; d<ndof; d++)
        {
            body2d->getGlobalIndex(i, d, globalIndex);
            ASSERT_EQ(counter, globalIndex);
            MatStencil stencil = {0, 0, i, d};
            body2d->getGlobalIndex(stencil, globalIndex);
            ASSERT_EQ(counter, globalIndex);
            counter++;
        }
}


TEST_F(SingleBodyTest, getGlobalIndex3D)
{
    PetscInt globalIndex, counter = 0;
    PetscInt ndof = 3;
    for (int i=0; i<body3d->nPts; i++)
        for (int d=0; d<ndof; d++)
        {
            body3d->getGlobalIndex(i, d, globalIndex);
            ASSERT_EQ(counter, globalIndex);
            MatStencil stencil = {0, 0, i, d};
            body3d->getGlobalIndex(stencil, globalIndex);
            ASSERT_EQ(counter, globalIndex);
            counter++;
        }
}


TEST_F(SingleBodyTest, calculateAvgForces2D)
{
    Vec f;
    DMCreateGlobalVector(body2d->da, &f);
    VecSet(f, 1.0);
    type::RealVec1D avg(2);
    body2d->calculateAvgForces(f, avg);
    ASSERT_EQ(body2d->dim, (int)avg.size());
    for (unsigned int d=0; d<avg.size(); d++)
        ASSERT_EQ(body2d->nPts * -1.0, avg[d]);
    VecDestroy(&f);
}


TEST_F(SingleBodyTest, calculateAvgForces3D)
{
    Vec f;
    DMCreateGlobalVector(body3d->da, &f);
    VecSet(f, 1.0);
    type::RealVec1D avg(3);
    body3d->calculateAvgForces(f, avg);
    ASSERT_EQ(body3d->dim, (int)avg.size());
    for (unsigned int d=0; d<avg.size(); d++)
        ASSERT_EQ(body3d->nPts * -1.0, avg[d]);
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
