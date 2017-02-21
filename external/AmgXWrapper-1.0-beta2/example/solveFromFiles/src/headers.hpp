/**
 * @file headers.hpp
 * @brief headers of Poisson3D test
 * @author Pi-Yueh Chuang (pychuang@gwu.edu)
 * @version alpha
 * @date 2015-11-06
 */
# pragma once

# define CHK CHKERRQ(ierr)

# include <iostream>
# include <algorithm>
# include <string>
# include <vector>
# include <cmath>
# include <cstring>

# include <petscsys.h>
# include <petsctime.h>
# include <petscksp.h>

# include "StructArgs.hpp"
# include "AmgXSolver.hpp"


PetscErrorCode readVec(Vec &vec, char *FN, const char *name);

PetscErrorCode readMat(Mat &mat, char *FN, const char *name);

PetscErrorCode createKSP(KSP &ksp, Mat &A, char *FN);

PetscErrorCode solve(KSP &ksp, Mat &A, Vec &u, Vec &rhs, Vec &u_exact, Vec &err,
        StructArgs &args, PetscLogEvent &warmUpEvent, PetscLogEvent &solvingEvent);

PetscErrorCode solve(AmgXSolver &amgx, Mat &A, Vec &u, Vec &rhs, Vec &u_exact, Vec &err,
        StructArgs &args, PetscLogEvent &warmUpEvent, PetscLogEvent &solvingEvent);

