/***************************************************************************//**
 * \file CartesianMeshTest.cpp
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Unit-test for the generation of uniform and stetched grids.
 */


#include "CartesianMesh.h"
#include "gtest/gtest.h"


class CartesianMeshTest : public ::testing::Test
{
public:
  CartesianMesh uniformMesh2d,
                uniformMesh3d,
                stretchedMesh3d;

  CartesianMeshTest()
  {
    std::string filePath = "CartesianMesh/cases/uniformMesh2d/cartesianMesh.yaml";
    uniformMesh2d = CartesianMesh(filePath);
    filePath = "CartesianMesh/cases/uniformMesh3d/cartesianMesh.yaml";
    uniformMesh3d = CartesianMesh(filePath);
    filePath = "CartesianMesh/cases/stretchedMesh3d/cartesianMesh.yaml";
    stretchedMesh3d = CartesianMesh(filePath);
  }
};

TEST_F(CartesianMeshTest, uniformMesh2dSize)
{
  EXPECT_EQ(uniformMesh2d.nx, 5);
  EXPECT_EQ(uniformMesh2d.ny, 4);
  EXPECT_EQ(uniformMesh2d.nz, 0);
}

TEST_F(CartesianMeshTest, uniformMesh2dCoordinates)
{
  PetscReal start;

  start = 0.0;
  for (auto i=uniformMesh2d.x.begin(); i!=uniformMesh2d.x.end(); i++)
  {
    EXPECT_DOUBLE_EQ(*i, start);
    start+=1.0;
  }
  
  for (auto i=uniformMesh2d.dx.begin(); i!=uniformMesh2d.dx.end(); i++)
  {
    EXPECT_DOUBLE_EQ(*i, 1);
  }

  start = 0.0;
  for (auto i=uniformMesh2d.y.begin(); i!=uniformMesh2d.y.end(); i++)
  {
    EXPECT_DOUBLE_EQ(*i, start);
    start+=1.0;
  }
  
  for (auto i=uniformMesh2d.dy.begin(); i!=uniformMesh2d.dy.end(); i++)
  {
    EXPECT_DOUBLE_EQ(*i, 1);
  }
}

TEST_F(CartesianMeshTest, uniformMesh3dSize)
{
  EXPECT_EQ(uniformMesh3d.nx, 5);
  EXPECT_EQ(uniformMesh3d.ny, 4);
  EXPECT_EQ(uniformMesh3d.nz, 6);
}

TEST_F(CartesianMeshTest, uniformMesh3dCoordinates)
{
  PetscReal start;

  start = 0.0;
  for (auto i=uniformMesh3d.x.begin(); i!=uniformMesh3d.x.end(); i++)
  {
    EXPECT_DOUBLE_EQ(*i, start);
    start+=1.0;
  }
  
  for (auto i=uniformMesh3d.dx.begin(); i!=uniformMesh3d.dx.end(); i++)
  {
    EXPECT_DOUBLE_EQ(*i, 1);
  }

  start = 0.0;
  for (auto i=uniformMesh3d.y.begin(); i!=uniformMesh3d.y.end(); i++)
  {
    EXPECT_DOUBLE_EQ(*i, start);
    start+=1.0;
  }
  
  for (auto i=uniformMesh3d.dy.begin(); i!=uniformMesh3d.dy.end(); i++)
  {
    EXPECT_DOUBLE_EQ(*i, 1);
  }

  start = 0.0;
  for (auto i=uniformMesh3d.z.begin(); i!=uniformMesh3d.z.end(); i++)
  {
    EXPECT_DOUBLE_EQ(*i, start);
    start+=1.0;
  }
  
  for (auto i=uniformMesh3d.dz.begin(); i!=uniformMesh3d.dz.end(); i++)
  {
    EXPECT_DOUBLE_EQ(*i, 1);
  }
}

TEST_F(CartesianMeshTest, stretchedMesh3dSize)
{
  EXPECT_EQ(stretchedMesh3d.nx, 3);
  EXPECT_EQ(stretchedMesh3d.ny, 4);
  EXPECT_EQ(stretchedMesh3d.nz, 5);
}

TEST_F(CartesianMeshTest, stretchedMesh3dCellWidths)
{
  PetscReal cellWidth;

  cellWidth = 1.0;
  for (auto i=stretchedMesh3d.dx.begin(); i!=stretchedMesh3d.dx.end(); i++)
  {
    EXPECT_DOUBLE_EQ(*i, cellWidth) << "The value should be " << cellWidth << " but is " << *i << std::endl;
    cellWidth*=1.5;
  }

  cellWidth = 1.0;
  for (auto i=stretchedMesh3d.dy.begin(); i!=stretchedMesh3d.dy.end(); i++)
  {
    EXPECT_DOUBLE_EQ(*i, cellWidth) << "The value should be " << cellWidth << " but is " << *i << std::endl;
    cellWidth*=2.0;
  }

  cellWidth = 1.0;
  for (auto i=stretchedMesh3d.dz.begin(); i!=stretchedMesh3d.dz.end(); i++)
  {
    EXPECT_DOUBLE_EQ(*i, cellWidth) << "The value should be " << cellWidth << " but is " << *i << std::endl;
    cellWidth*=1.1;
  }
}

int main(int argc, char **argv)
{
  PetscErrorCode ierr, result;

  ::testing::InitGoogleTest(&argc, argv);
  ierr = PetscInitialize(&argc, &argv, NULL, NULL); CHKERRQ(ierr);
  result = RUN_ALL_TESTS();
  ierr = PetscFinalize(); CHKERRQ(ierr);

  return result;
}