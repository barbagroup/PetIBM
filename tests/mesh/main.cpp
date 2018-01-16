/*
 * main.cpp
 * Copyright (C) 2017 Pi-Yueh Chuang <pychuang@gwu.edu>
 *
 * Distributed under terms of the MIT license.
 */

# include <petsc.h>
# include <gtest/gtest.h>

# include <iostream>
# include <sstream>

class PETScEnvironment : public ::testing::Environment
{
    public:
        
        virtual ~PETScEnvironment() {};
        virtual void SetUp()
        {
            char **argv;
            int argc=0;
            PetscErrorCode ierr = PetscInitialize(&argc, &argv, nullptr, nullptr);
            ASSERT_FALSE(ierr);
        }
        
        virtual void TearDown()
        {
            PetscErrorCode ierr = PetscFinalize();
            ASSERT_FALSE(ierr);
        }
};


// Run all tests
int main(int argc, char **argv)
{
    PetscErrorCode status;

    ::testing::InitGoogleTest(&argc, argv);
    ::testing::AddGlobalTestEnvironment(new PETScEnvironment);
    
    status = RUN_ALL_TESTS();

    return status;
} // main
