/**
 * \file fixSingularMat.hpp
 * \brief prototypes of fixSingularMat.
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \date 2017-06-01
 */


# pragma once

// PETSc
# include <petscsys.h>
# include <petscmat.h>
# include <petscvec.h>


/**
 * \brief modify matrix A and RHS vector for all-Neumann BCs to avoid sungular matrix.
 *
 * All-Neumann BC causes issues like indefinite matrix or infinite number of 
 * solutions. We hence assume the exact solution at the node i = j = 0 is known.
 * (The same way we deal with the pressure in CFD applications.) The entry 
 * A[0, 0] is modified to a value close to other diagonal coefficients, and
 * A[0, *] = A[*, 0] = 0; also all corresponding entries in RHS vector is
 * modified respectively.
 *
 * In other words, we apply a Dirichlet BC at one point.
 *
 * \param A [in, out] matrix A generated from the function generateA(...)
 * \param RHS [in, out] the right-hand-side vector generated from generateRHS(...)
 * \param exact [in] the exact solution vector generated from generateExt(...)
 *
 * \return PetscErrorCode.
 */
PetscErrorCode fixSingularMat(Mat &A, Vec &RHS, const Vec &exact);
