/**
 * \file delta_test.cpp
 * \brief Runs unit-tests for the PetIBM utilities.
 */

#include <petsc.h>

#include "gtest/gtest.h"

#include "delta_test.inl"
#include "singlebody_test.inl"


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
