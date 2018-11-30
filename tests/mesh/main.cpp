/*
 * \file tests/mesh/main.cpp
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 */

#include <gtest/gtest.h>
#include <petsc.h>

// Run all tests
int main(int argc, char **argv)
{
    PetscErrorCode status;

    ::testing::InitGoogleTest(&argc, argv);
    status = PetscInitialize(&argc, &argv, nullptr, nullptr);
    CHKERRQ(status);
    status = RUN_ALL_TESTS();
    status = PetscFinalize();
    CHKERRQ(status);

    return status;
}  // main
