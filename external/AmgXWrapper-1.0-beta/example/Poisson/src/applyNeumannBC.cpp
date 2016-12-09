/**
 * @file applyNeumannBC.cpp
 * @brief Definition of function applyNeumannBC(...)
 * @author Pi-Yueh Chuang (pychuang@gwu.edu)
 * @version beta
 * @date 2016-02-13
 */
# include "headers.hpp"

/**
 * @brief Modify matrix A and RHS vector for all-Neumann BCs
 *
 * All Neumann BCs causes issues like indefinite matrix or infinite number of 
 * solutions. We hence assume the exact solution at the node i = j = 0 is known.
 * (The same way we deal with the pressure in CFD applications.) The entry 
 * A[0, 0] is modified to 1, and A[0, *] = A[*, 0] = 0; also all corresponding
 * entries in RHS vector is modified respectively.
 *
 * @param A Matrix A generated from the function generateA(...)
 * @param RHS The right-hand-side vector generated from generateRHS(...)
 * @param exact The exact solution vector generated from generateExt(...)
 *
 * @return 
 */
PetscErrorCode applyNeumannBC(Mat &A, Vec &RHS, const Vec &exact)
{
    PetscErrorCode      ierr;

    PetscInt            row[1] = {0};

    ierr = MatZeroRowsColumns(A, 1, row, 1.0, exact, RHS);         CHKERRQ(ierr);

    return 0;
}
