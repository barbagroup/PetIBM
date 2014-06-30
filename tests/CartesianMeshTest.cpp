#include "CartesianMesh.h"
#include "gtest/gtest.h"

class CartesianMeshTest : public ::testing::Test
{
protected:

	CartesianMesh UniformMesh2d,
	              UniformMesh3d,
	              StretchedMesh3d;

	CartesianMeshTest()
	{
		UniformMesh2d.initialize("tests/CartesianMesh/UniformMesh2d.yaml");
		UniformMesh3d.initialize("tests/CartesianMesh/UniformMesh3d.yaml");
		StretchedMesh3d.initialize("tests/CartesianMesh/StretchedMesh3d.yaml");
	}
};

TEST_F(CartesianMeshTest, UniformMesh2dSize)
{
	EXPECT_EQ(UniformMesh2d.nx, 5);
	EXPECT_EQ(UniformMesh2d.ny, 4);
	EXPECT_EQ(UniformMesh2d.nz, 0);
}

TEST_F(CartesianMeshTest, UniformMesh2dCoordinates)
{
	PetscReal start;

	start = 0.0;
	for(auto i=UniformMesh2d.x.begin(); i!=UniformMesh2d.x.end(); i++)
	{
		EXPECT_DOUBLE_EQ(*i, start);
		start+=1.0;
	}
	
	for(auto i=UniformMesh2d.dx.begin(); i!=UniformMesh2d.dx.end(); i++)
	{
		EXPECT_DOUBLE_EQ(*i, 1);
	}

	start = 0.0;
	for(auto i=UniformMesh2d.y.begin(); i!=UniformMesh2d.y.end(); i++)
	{
		EXPECT_DOUBLE_EQ(*i, start);
		start+=1.0;
	}
	
	for(auto i=UniformMesh2d.dy.begin(); i!=UniformMesh2d.dy.end(); i++)
	{
		EXPECT_DOUBLE_EQ(*i, 1);
	}
}

TEST_F(CartesianMeshTest, UniformMesh3dSize)
{
	EXPECT_EQ(UniformMesh3d.nx, 5);
	EXPECT_EQ(UniformMesh3d.ny, 4);
	EXPECT_EQ(UniformMesh3d.nz, 6);
}

TEST_F(CartesianMeshTest, UniformMesh3dCoordinates)
{
	PetscReal start;

	start = 0.0;
	for(auto i=UniformMesh3d.x.begin(); i!=UniformMesh3d.x.end(); i++)
	{
		EXPECT_DOUBLE_EQ(*i, start);
		start+=1.0;
	}
	
	for(auto i=UniformMesh3d.dx.begin(); i!=UniformMesh3d.dx.end(); i++)
	{
		EXPECT_DOUBLE_EQ(*i, 1);
	}

	start = 0.0;
	for(auto i=UniformMesh3d.y.begin(); i!=UniformMesh3d.y.end(); i++)
	{
		EXPECT_DOUBLE_EQ(*i, start);
		start+=1.0;
	}
	
	for(auto i=UniformMesh3d.dy.begin(); i!=UniformMesh3d.dy.end(); i++)
	{
		EXPECT_DOUBLE_EQ(*i, 1);
	}

	start = 0.0;
	for(auto i=UniformMesh3d.z.begin(); i!=UniformMesh3d.z.end(); i++)
	{
		EXPECT_DOUBLE_EQ(*i, start);
		start+=1.0;
	}
	
	for(auto i=UniformMesh3d.dz.begin(); i!=UniformMesh3d.dz.end(); i++)
	{
		EXPECT_DOUBLE_EQ(*i, 1);
	}
}

TEST_F(CartesianMeshTest, StretchedMesh3dSize)
{
	EXPECT_EQ(StretchedMesh3d.nx, 3);
	EXPECT_EQ(StretchedMesh3d.ny, 4);
	EXPECT_EQ(StretchedMesh3d.nz, 5);
}

TEST_F(CartesianMeshTest, StretchedMesh3dCellWidths)
{
	PetscReal cellWidth;

	cellWidth = 1.0;
	for(auto i=StretchedMesh3d.dx.begin(); i!=StretchedMesh3d.dx.end(); i++)
	{
		EXPECT_DOUBLE_EQ(*i, cellWidth) << "The value should be " << cellWidth << " but is " << *i << std::endl;
		cellWidth*=1.5;
	}

	cellWidth = 1.0;
	for(auto i=StretchedMesh3d.dy.begin(); i!=StretchedMesh3d.dy.end(); i++)
	{
		EXPECT_DOUBLE_EQ(*i, cellWidth) << "The value should be " << cellWidth << " but is " << *i << std::endl;
		cellWidth*=2.0;
	}

	cellWidth = 1.0;
	for(auto i=StretchedMesh3d.dz.begin(); i!=StretchedMesh3d.dz.end(); i++)
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
