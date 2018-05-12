/**
 * \file createKSP.hpp
 * \brief prototype of createKSP
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \date 2017-05-31
 */


# pragma once


// PETSc
# include <petscsys.h>
# include <petscdmda.h>
# include <petscmat.h>
# include <petscvec.h>
# include <petscksp.h>


/**
 * \brief create a KSP solver based on settings in a file.
 *
 * \param ksp [out] the resulting KSP solver.
 * \param A [in] the coefficient bound to the solver.
 * \param grid [in] a PETSc DMDA to be bound wiht KSP.
 * \param FN [in] a ASCII file containing KSP settings.
 *
 * \return PetscErrorCode.
 */
PetscErrorCode createKSP(KSP &ksp, Mat &A, DM &grid, char *FN);
