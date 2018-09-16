/**
 * \file solution.cpp
 * \brief Implementations of createSolution and members of SolutionBase.
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
}  // ~SolutionBase

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
}  // destroy

PetscErrorCode SolutionBase::printInfo() const
{
    PetscFunctionBeginUser;

    PetscErrorCode ierr;

    ierr = io::print(info); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // printInfo

PetscErrorCode createSolution(const type::Mesh &mesh, type::Solution &solution)
{
    PetscFunctionBeginUser;

    solution = std::make_shared<SolutionSimple>(mesh);

    PetscFunctionReturn(0);
}  // createSolution

}  // end of namespace solution
}  // end of namespace petibm
