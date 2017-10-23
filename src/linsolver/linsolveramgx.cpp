/*! Implementation of the methods of the class `LinSolverAmgX`.
 * \file amgxsolver.cpp
 */


// PetIBM
#include <petibm/linsolveramgx.h>


namespace petibm
{
namespace linsolver
{

// constructor
LinSolverAmgX::LinSolverAmgX(const std::string &_name, const std::string &_config):
    LinSolverBase(_name, _config) { init(); }


// destructor
LinSolverAmgX::~LinSolverAmgX() { amgx.finalize(); }


// underlying initialization function
PetscErrorCode LinSolverAmgX::init()
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;
    
    type = "NVIDIA AmgX";

    ierr = amgx.initialize(PETSC_COMM_WORLD, "dDDI", config); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


// set coefficient matrix
PetscErrorCode LinSolverAmgX::setMatrix(const Mat &A)
{
    PetscErrorCode ierr;

    ierr = amgx.setA(A); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // setMatrix


// solve linear system
PetscErrorCode LinSolverAmgX::solve(Vec &x, Vec &b)
{
    PetscErrorCode ierr;

    ierr = amgx.solve(x, b); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // solve


// get number of iterations
PetscErrorCode LinSolverAmgX::getIters(PetscInt &iters)
{
    PetscErrorCode ierr;

    iters = amgx.getIters();

    PetscFunctionReturn(0);
} // getIters

} // end of namespace linsolver
} // end of namespace petibm
