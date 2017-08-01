/**
 * \file delta_test.cpp
 * \brief Unit-tests for the discrete delta functions.
 */

#include <vector>
#include <random>
#include <algorithm>

#include "gtest/gtest.h"

#include "utilities/delta.h"

using namespace petibm::utilities::delta;


// check value is zero outside region of influence
TEST(utilitiesDeltaRomaEtAlTest, zeroOutside)
{
	const PetscReal h = 1.0;
	EXPECT_EQ(0.0, Roma_et_al(1.5, h));
	EXPECT_EQ(0.0, Roma_et_al(2.0, h));
}

// check maximum value at 0
TEST(utilitiesDeltaRomaEtAlTest, maximumValue)
{
	const PetscReal h = 1.0;
	EXPECT_EQ(2.0/3.0, Roma_et_al(0.0, h));
}

// check delta function is monotonically decreasing
TEST(utilitiesDeltaRomaEtAlTest, decreasingInfluence)
{
	const PetscReal h = 1.0;
	// create a sorted vector of random reals between 0.0 and 1.5
	std::vector<PetscReal> vals(10);
	std::default_random_engine engine;
	std::uniform_real_distribution<PetscReal> distrib(0.0, 1.5);
	std::generate(vals.begin(), vals.end(), [&](){ return distrib(engine); });
	std::sort(vals.begin(), vals.end());
	// assert decreasing influence as the distance increases
	for (unsigned int i=0; i<vals.size()-1; i++)
		ASSERT_GT(Roma_et_al(vals[i], h), Roma_et_al(vals[i+1], h));
}
