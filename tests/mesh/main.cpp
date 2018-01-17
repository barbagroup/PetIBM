/*
 * main.cpp
 * Copyright (C) 2017 Pi-Yueh Chuang <pychuang@gwu.edu>
 *
 * Distributed under terms of the MIT license.
 */

# include <petsc.h>
# include <gtest/gtest.h>


// Run all tests
int main(int argc, char **argv)
{
    PetscErrorCode status;

    ::testing::InitGoogleTest(&argc, argv);
    status = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(status);
    status = RUN_ALL_TESTS();
    status = PetscFinalize(); CHKERRQ(status);

    return status;
} // main
