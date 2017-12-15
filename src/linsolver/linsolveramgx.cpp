/** 
 * \file linsolveramgx.cpp
 * \brief Implementation of LinSolverAmgX.
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \copyright MIT.
 */


// PetIBM
#include <petibm/linsolveramgx.h>


std::vector<AmgXSolver*> solvers;

PetscErrorCode preDestroyAmgXSolvers()
{
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    for(auto &it: solvers)
    {
        ierr = it->finalize(); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
}

namespace petibm
{
namespace linsolver
{

// constructor
LinSolverAmgX::LinSolverAmgX(const std::string &_name, const std::string &_config):
    LinSolverBase(_name, _config) { init(); }


// destructor
LinSolverAmgX::~LinSolverAmgX()
{ 
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    PetscBool finalized;

    ierr = PetscFinalized(&finalized); CHKERRV(ierr);
    if (finalized) return;

    ierr = amgx.finalize(); CHKERRV(ierr);
}


// manually destroy
PetscErrorCode LinSolverAmgX::destroy()
{
    PetscFunctionBeginUser;

    PetscErrorCode ierr;

    ierr = amgx.finalize(); CHKERRQ(ierr);
    ierr = LinSolverBase::destroy(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


// underlying initialization function
PetscErrorCode LinSolverAmgX::init()
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;
    
    type = "NVIDIA AmgX";

    ierr = amgx.initialize(PETSC_COMM_WORLD, "dDDI", config); CHKERRQ(ierr);

    solvers.push_back(&amgx);
    ierr = PetscRegisterFinalize(&preDestroyAmgXSolvers); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


// set coefficient matrix
PetscErrorCode LinSolverAmgX::setMatrix(const Mat &A)
{
    PetscFunctionBeginUser;

    PetscErrorCode ierr;

    ierr = amgx.setA(A); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // setMatrix


// solve linear system
PetscErrorCode LinSolverAmgX::solve(Vec &x, Vec &b)
{
    PetscFunctionBeginUser;

    PetscErrorCode ierr;

    ierr = amgx.solve(x, b); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // solve


// get number of iterations
PetscErrorCode LinSolverAmgX::getIters(PetscInt &iters)
{
    PetscFunctionBeginUser;

    PetscErrorCode ierr;

    ierr = amgx.getIters(iters); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // getIters


// get the norm of final residual
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
