/**
 * \file solution.cpp
 * \brief Implementation of the factory function
 *        and members of the class petibm::solution::SolutionBase.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 */

#include <petibm/io.h>
#include <petibm/solution.h>
#include <petibm/solutionsimple.h>

namespace petibm
{
namespace solution
{
// Destructor.
SolutionBase::~SolutionBase()
{
    PetscErrorCode ierr;
    PetscBool finalized;

    PetscFunctionBeginUser;

    ierr = PetscFinalized(&finalized); CHKERRV(ierr);
    if (finalized) return;

    ierr = VecDestroy(&UGlobal); CHKERRV(ierr);
    ierr = VecDestroy(&pGlobal); CHKERRV(ierr);
    comm = MPI_COMM_NULL;
}  // ~SolutionBase

//  Manually destroy data.
PetscErrorCode SolutionBase::destroy()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    dim = -1;
    ierr = VecDestroy(&UGlobal); CHKERRQ(ierr);
    ierr = VecDestroy(&pGlobal); CHKERRQ(ierr);
    info = "";

    comm = MPI_COMM_NULL;
    mpiRank = mpiSize = 0;
    mesh.reset();

    PetscFunctionReturn(0);
}  // destroy

// Print information about the solution to standard output.
PetscErrorCode SolutionBase::printInfo() const
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = io::print(info); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // printInfo

// Factory function to create a petibm::solution::Solution object.
PetscErrorCode createSolution(const type::Mesh &mesh, type::Solution &solution)
{
    PetscFunctionBeginUser;

    solution = std::make_shared<SolutionSimple>(mesh);

    PetscFunctionReturn(0);
}  // createSolution

}  // end of namespace solution
}  // end of namespace petibm
