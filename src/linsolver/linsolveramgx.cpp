/*! Implementation of the methods of the class `LinSolverAmgX`.
 * \file amgxsolver.cpp
 */


#include <iostream>
#include <fstream>
#include <cstdio>

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

    // create temporary empty file if no configuration file is provided
    char fname[L_tmpnam];
    if (config == "None")
    {
    	config = std::tmpnam(fname);
    	std::fstream file(fname, std::ios::out);
    }

  	ierr = PetscPrintf(PETSC_COMM_WORLD,
                       "[%s solver] Configuration file: %s\n",
                       name.c_str(), config.c_str()); CHKERRQ(ierr);

    ierr = amgx.initialize(PETSC_COMM_WORLD, "dDDI", config); CHKERRQ(ierr);

    solvers.push_back(&amgx);
    ierr = PetscRegisterFinalize(&preDestroyAmgXSolvers); CHKERRQ(ierr);

    // remove temporary file if necessary
    if (config == fname)
    	std::remove(config.c_str());

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
