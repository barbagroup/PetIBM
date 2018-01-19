/** 
 * \file linsolveramgx.cpp
 * \brief Implementation of LinSolverAmgX.
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \copyright MIT.
 */


// PetIBM
#include <petibm/linsolveramgx.h>


// a global (in this file scope) vector of references to AmgXSolvers
std::vector<AmgXSolver*> solvers;

// a wrapper function for PetscFinalize() to destroy AmgXSolvers automatically
PetscErrorCode preDestroyAmgXSolvers()
{
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    for(auto &it: solvers)
    {
        ierr = it->finalize(); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
} // preDestroyAmgXSolvers

namespace petibm
{
namespace linsolver
{

// implement LinSolverAmgX::LinSolverAmgX
LinSolverAmgX::LinSolverAmgX(
        const std::string &solverName, const std::string &file):
    LinSolverBase(solverName, file)
{
	init();
} // LinSolverAmgX


// implement LinSolverAmgX::~LinSolverAmgX
LinSolverAmgX::~LinSolverAmgX()
{ 
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    PetscBool finalized;

    ierr = PetscFinalized(&finalized); CHKERRV(ierr);
    if (finalized) return;

    ierr = amgx.finalize(); CHKERRV(ierr);
} // ~LinSolverAmgX


// implement LinSolverAmgX::destroy
PetscErrorCode LinSolverAmgX::destroy()
{
    PetscFunctionBeginUser;

    PetscErrorCode ierr;

    ierr = amgx.finalize(); CHKERRQ(ierr);
    ierr = LinSolverBase::destroy(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // destroy


// implement LinSolverAmgX::init
PetscErrorCode LinSolverAmgX::init()
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;
    
    type = "NVIDIA AmgX";

    ierr = amgx.initialize(PETSC_COMM_WORLD, "dDDI", config); CHKERRQ(ierr);

    solvers.push_back(&amgx);
    ierr = PetscRegisterFinalize(&preDestroyAmgXSolvers); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // init


// implement LinSolverAmgX::setMatrix
PetscErrorCode LinSolverAmgX::setMatrix(const Mat &A)
{
    PetscFunctionBeginUser;

    PetscErrorCode ierr;

    ierr = amgx.setA(A); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // setMatrix


// implement LinSolverAmgX::solve
PetscErrorCode LinSolverAmgX::solve(Vec &x, Vec &b)
{
    PetscFunctionBeginUser;

    PetscErrorCode ierr;

    ierr = amgx.solve(x, b); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // solve


// implement LinSolverAmgX::getIters
PetscErrorCode LinSolverAmgX::getIters(PetscInt &iters)
{
    PetscFunctionBeginUser;

    PetscErrorCode ierr;

    ierr = amgx.getIters(iters); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // getIters


// implement LinSolverAmgX::getResidual
PetscErrorCode LinSolverAmgX::getResidual(PetscReal &res)
{
    PetscFunctionBeginUser;

    PetscErrorCode ierr;
    
    PetscInt iter;
    
    ierr = amgx.getIters(iter); CHKERRQ(ierr);
    ierr = amgx.getResidual(iter-1, res); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // getResidual

} // end of namespace linsolver
} // end of namespace petibm
