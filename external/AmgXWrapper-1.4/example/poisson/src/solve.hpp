/**
 * \file solve.hpp
 * \brief prototypes of solve.
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \date 2017-05-31
 */


# pragma once


// PETSc
# include <petscsys.h>
# include <petscmat.h>
# include <petscvec.h>
# include <petscksp.h>

// AmgXWrapper
# include "AmgXSolver.hpp"

// header
# include "StructArgs.hpp"


/**
 * \brief solve the system with KSP and return performance info.
 *
 * \param ksp [in] KSP solver.
 * \param A [in] cosfficient matrix.
 * \param lhs [out] left hand side.
 * \param rhs [in] right hand side.
 * \param exact [in] exact solution for error check.
 * \param err [out] errors.
 * \param args [in] a StructArgs instance.
 * \param warmUpEvent [out] PetscLogEvent for warming up cycle.
 * \param solvingEvent [out] PetscLogEvent for solving.
 *
 * \return PetscErrorCode.
 */
PetscErrorCode solve(KSP &ksp, Mat &A, Vec &lhs, Vec &rhs, Vec &exact, Vec &err,
        StructArgs &args, PetscLogEvent &warmUpEvent, PetscLogEvent &solvingEvent);


/**
 * \brief solve the system with KSP and return performance info.
 *
 * \param amgx [in] AmgXSolver instance.
 * \param A [in] cosfficient matrix.
 * \param lhs [out] left hand side.
 * \param rhs [in] right hand side.
 * \param exact [in] exact solution for error check.
 * \param err [out] errors.
 * \param args [in] a StructArgs instance.
 * \param warmUpEvent [out] PetscLogEvent for warming up cycle.
 * \param solvingEvent [out] PetscLogEvent for solving.
 *
 * \return PetscErrorCode.
 */
PetscErrorCode solve(AmgXSolver &amgx, Mat &A, Vec &lhs, Vec &rhs, Vec &exact, Vec &err,
        StructArgs &args, PetscLogEvent &warmUpEvent, PetscLogEvent &solvingEvent);
