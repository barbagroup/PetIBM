/*
 * solution.cpp
 * Copyright (C) 2017 Pi-Yueh Chuang <pychuang@gwu.edu>
 *
 * Distributed under terms of the MIT license.
 */

# include <petibm/solution.h>
# include <petibm/solutionsimple.h>
# include <petibm/io.h>


namespace petibm
{
namespace solution
{
    PetscErrorCode SolutionBase::printInfo() const
    {
        PetscFunctionBeginUser;
        
        PetscErrorCode  ierr;
        
        ierr = io::print(info); CHKERRQ(ierr);
        
        PetscFunctionReturn(0);
    }
    

    PetscErrorCode createSolution(
            const type::Mesh &mesh, type::Solution &solution)
    {
        PetscFunctionBeginUser;
        
        solution = std::make_shared<SolutionSimple>(mesh);

        PetscFunctionReturn(0);
    }
}
}
