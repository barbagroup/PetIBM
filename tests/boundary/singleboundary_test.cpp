/**
 * \file singleboundary_test.cpp
 * \brief Unit-tests for the class `SingleBoundary`.
 */

#include <petsc.h>

#include <gtest/gtest.h>
#include <yaml-cpp/yaml.h>

#include <petibm/singleboundary.h>
#include <petibm/parser.h>
#include <petibm/mesh.h>

using namespace petibm;


class SingleBoundaryTest : public ::testing::Test
{
protected:

    SingleBoundaryTest(){  };

    virtual ~SingleBoundaryTest(){  };

    virtual void SetUp()
    {
        using namespace YAML;
        
        Node config;
        type::Mesh  mesh;
        
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
        petibm::boundary::createSingleBoundary(mesh, type::BCLoc::XMINUS, 
                type::Field::u, 0.0, type::BCType::DIRICHLET, boundary);
    };

    virtual void TearDown(){  };

    type::SingleBoundary boundary;
}; // SingleBoundaryTest


TEST_F(SingleBoundaryTest, init)
{
    ASSERT_EQ(3, boundary->dim);
    ASSERT_EQ(type::BCLoc::XMINUS, boundary->loc);
    PetscMPIInt size;
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    if (size == 1)
        ASSERT_EQ(PETSC_TRUE, boundary->onThisProc);
    ASSERT_EQ(type::BCType::DIRICHLET, boundary->type);
    ASSERT_EQ(0.0, boundary->value);
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
