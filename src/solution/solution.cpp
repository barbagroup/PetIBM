/*
 * solution.cpp
 * Copyright (C) 2017 Pi-Yueh Chuang <pychuang@gwu.edu>
 *
 * Distributed under terms of the MIT license.
 */

# include <petibm/solution.h>
# include <petibm/solutionsimple.h>


namespace petibm
{
namespace solution
{
    PetscErrorCode createSolution(
            const type::Mesh &mesh, type::Solution &solution)
    {
        PetscFunctionBeginUser;
        
        solution = std::make_shared<SolutionSimple>(mesh);

        PetscFunctionReturn(0);
    }
}
}
