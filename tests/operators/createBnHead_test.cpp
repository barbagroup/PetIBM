/**
 * \file createBnHead_test.cpp
 * \brief Unit-tests for the function `petibm::operators::createBnHead`.
 */

#include <cmath>

#include <petsc.h>

#include "gtest/gtest.h"

#include "operators/operators.h"


// Test Mat BNHat up to N=10 when Op is a diagonal Mat
TEST(operatorsCreateBnHeadTest, testBNHatDiagonalOp)
{
	Mat Op;  // operator (e.g., Laplacian)
	PetscReal dt = 2.3,  // time-step size
	          c = 0.5,  // time-scheme coefficient of the implicit diffusive term
	          val = 2.0/dt; // value set on the diagonal of the operator
	PetscInt nx = 10,  // number of points in the x-direction
	         ny = 12;  // number of points in the y-direction
	PetscReal ans = nx*ny*dt;  // expected sum of all elements of B1Hat
	// Create and assemble operator
	MatCreate(PETSC_COMM_WORLD, &Op);
	MatSetSizes(Op, PETSC_DECIDE, PETSC_DECIDE, nx*ny, nx*ny);
	MatSetType(Op, MATAIJ);
	MatSetUp(Op);
	for (PetscInt i=0; i<nx*ny; i++)
		MatSetValues(Op, 1, &i, 1, &i, &val, INSERT_VALUES);
	MatAssemblyBegin(Op, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(Op, MAT_FINAL_ASSEMBLY);
	for (PetscInt N=1; N<=10; N++)
	{
		Mat BNHat;  // Nth-order approximation of the inverse of (I/dt - c*Op)
		// Call function to test
		petibm::operators::createBnHead(Op, dt, c, N, BNHat);
		// Check size of Mat BNHat
		{
			PetscInt nrows, ncols;
			MatGetSize(BNHat, &nrows, &ncols);
			ASSERT_EQ(nx*ny, nrows);
			ASSERT_EQ(nx*ny, ncols);
		}
		// Check sum of elements of BNHat is the expected value
		if (N > 1)
			ans += dt * nx*ny * std::pow(c*dt*val, N-1);
		{
			PetscReal sum;
			Vec v;
			MatCreateVecs(Op, &v, nullptr);
			MatGetRowSum(BNHat, v);
			VecSum(v, &sum);
			ASSERT_TRUE(abs(sum - ans) <= 1.0E-16);
			VecDestroy(&v);
		}
		MatDestroy(&BNHat);
	}
	MatDestroy(&Op);
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
