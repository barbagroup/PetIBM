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
    SolutionBase::~SolutionBase()
    {
        PetscFunctionBeginUser;
        PetscErrorCode ierr;
        PetscBool finalized;

        ierr = PetscFinalized(&finalized); CHKERRV(ierr);
        if (finalized) return;

        ierr = VecDestroy(&UGlobal); CHKERRV(ierr);
        ierr = VecDestroy(&pGlobal); CHKERRV(ierr);
        comm = MPI_COMM_NULL;
    }


    PetscErrorCode SolutionBase::destroy()
    {
        PetscFunctionBeginUser;
        PetscErrorCode ierr;

        dim = -1;
        ierr = VecDestroy(&UGlobal); CHKERRQ(ierr);
        ierr = VecDestroy(&pGlobal); CHKERRQ(ierr);
        info = "";

        comm = MPI_COMM_NULL;
        mpiRank = mpiSize = 0;
        mesh.reset();

        PetscFunctionReturn(0);
    }


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
