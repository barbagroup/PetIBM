/**
 * \file cartesianmesh_test.cpp
 * \brief Unit-tests for the class `CartesianMesh`.
 */

#include <petsc.h>

#include <gtest/gtest.h>
#include <yaml-cpp/yaml.h>

#include <petibm/mesh.h>

using namespace petibm;


class CartesianMeshTest3D : public ::testing::Test
{
protected:

    CartesianMeshTest3D(){  };

    virtual ~CartesianMeshTest3D(){  };

    virtual void SetUp()
    {
        using namespace YAML;
        
        Node config;
        
        config["mesh"].push_back(Node(NodeType::Map));
        config["mesh"][0]["direction"] = "x";
        config["mesh"][1]["direction"] = "y";
        config["mesh"][2]["direction"] = "z";
        for(unsigned int i=0; i<3; ++i)
        {
            config["mesh"][i]["start"] = "0.1";
            config["mesh"][i]["subDomains"].push_back(Node(NodeType::Map));
            config["mesh"][i]["subDomains"][0]["end"] = 1.1;
            config["mesh"][i]["subDomains"][0]["cells"] = 10;
            config["mesh"][i]["subDomains"][0]["stretchRatio"] = 1.0;
        }
        
        config["flow"] = YAML::Node(NodeType::Map);
        config["flow"]["boundaryConditions"].push_back(Node(NodeType::Map));
        config["flow"]["boundaryConditions"][0]["location"] = "xMinus";
        config["flow"]["boundaryConditions"][1]["location"] = "xPlus";
        config["flow"]["boundaryConditions"][2]["location"] = "yMinus";
        config["flow"]["boundaryConditions"][3]["location"] = "yPlus";
        config["flow"]["boundaryConditions"][4]["location"] = "zMinus";
        config["flow"]["boundaryConditions"][5]["location"] = "zPlus";
        
        for(unsigned int i=0; i<4; ++i)
        {
            config["flow"]["boundaryConditions"][i]["u"][0] = "DIRICHLET";
            config["flow"]["boundaryConditions"][i]["u"][1] = 0.0;
            config["flow"]["boundaryConditions"][i]["v"][0] = "DIRICHLET";
            config["flow"]["boundaryConditions"][i]["v"][1] = 0.0;
            config["flow"]["boundaryConditions"][i]["w"][0] = "DIRICHLET";
            config["flow"]["boundaryConditions"][i]["w"][1] = 0.0;
        }
        config["flow"]["boundaryConditions"][3]["u"][1] = 1.0;
        
        for(unsigned int i=4; i<6; ++i)
        {
            config["flow"]["boundaryConditions"][i]["u"][0] = "PERIODIC";
            config["flow"]["boundaryConditions"][i]["u"][1] = 0.0;
            config["flow"]["boundaryConditions"][i]["v"][0] = "PERIODIC";
            config["flow"]["boundaryConditions"][i]["v"][1] = 0.0;
            config["flow"]["boundaryConditions"][i]["w"][0] = "PERIODIC";
            config["flow"]["boundaryConditions"][i]["w"][1] = 0.0;
        }
        
        petibm::mesh::createMesh(PETSC_COMM_WORLD, config, mesh);

    };

    virtual void TearDown(){  };

    type::Mesh mesh;
}; // CartesianMeshTest


TEST_F(CartesianMeshTest3D, init3D)
{
    // check dimension
    ASSERT_EQ(3, mesh->dim);
    // check domain limits
    for (unsigned int d=0; d<3; d++)
    {
        ASSERT_EQ(0.1, mesh->min[d]);
        ASSERT_EQ(1.1, mesh->max[d]);
    }
    // check number of points for each field
    std::vector<PetscReal> n_exp(3);
    // for x-velocity
    n_exp = {9, 10, 10};
    for (unsigned int d=0; d<3; d++)
        ASSERT_EQ(n_exp[d], mesh->n[0][d]);
    // for y-velocity
    n_exp = {10, 9, 10};
    for (unsigned int d=0; d<3; d++)
        ASSERT_EQ(n_exp[d], mesh->n[1][d]);
    // for z-velocity
    n_exp = {10, 10, 10};
    for (unsigned int d=0; d<3; d++)
        ASSERT_EQ(n_exp[d], mesh->n[2][d]);
    // for pressure
    n_exp = {10, 10, 10};
    for (unsigned int d=0; d<3; d++)
        ASSERT_EQ(n_exp[d], mesh->n[3][d]);
    // for vertices
    n_exp = {11, 11, 11};
    for (unsigned int d=0; d<3; d++)
        ASSERT_EQ(n_exp[d], mesh->n[4][d]);
}
