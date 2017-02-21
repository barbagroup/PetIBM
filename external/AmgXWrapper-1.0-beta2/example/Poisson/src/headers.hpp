/**
 * @file headers.hpp
 * @brief headers of Poisson3D test
 * @author Pi-Yueh Chuang (pychuang@gwu.edu)
 * @version alpha
 * @date 2015-11-06
 */
# pragma once

# include <iostream>
# include <algorithm>
# include <string>
# include <vector>
# include <cmath>
# include <cstring>

# include <petscsys.h>
# include <petsctime.h>
# include <petscksp.h>
# include <petscdmda.h>

# include "StructArgs.hpp"
# include "AmgXSolver.hpp"


# define CHK CHKERRQ(ierr)

# define c1 2.0*1.0*M_PI


PetscErrorCode generateGrid(const DM &grid, 
        const PetscInt &Nx, const PetscInt &Ny, const PetscInt &Nz,
        const PetscScalar &Lx, const PetscScalar &Ly, const PetscScalar &Lz,
        PetscScalar &dx, PetscScalar &dy, PetscScalar &dz,
        Vec &x, Vec &y, Vec &z);

PetscErrorCode generateGrid(const DM &grid, 
        const PetscInt &Nx, const PetscInt &Ny,
        const PetscScalar &Lx, const PetscScalar &Ly,
        PetscScalar &dx, PetscScalar &dy,
        Vec &x, Vec &y);

PetscErrorCode generateRHS(const DM &grid, 
        const Vec &x, const Vec &y, const Vec &z, Vec &rhs);

PetscErrorCode generateRHS(const DM &grid, 
        const Vec &x, const Vec &y, Vec &rhs);

PetscErrorCode generateExt(const DM &grid, 
        const Vec &x, const Vec &y, const Vec &z, Vec &exact);

PetscErrorCode generateExt(const DM &grid, 
        const Vec &x, const Vec &y, Vec &exact);

PetscErrorCode generateA(const DM &grid, 
        const PetscScalar &dx, const PetscScalar &dy, const PetscScalar &dz, Mat &A);

PetscErrorCode generateA(const DM &grid, 
        const PetscScalar &dx, const PetscScalar &dy, Mat &A);

PetscErrorCode createKSP(KSP &ksp, Mat &A, DM &grid, char *FN);

PetscErrorCode solve(KSP &ksp, Mat &A, Vec &u, Vec &rhs, Vec &u_exact, Vec &err,
        StructArgs &args, PetscLogEvent &warmUpEvent, PetscLogEvent &solvingEvent);

PetscErrorCode solve(AmgXSolver &amgx, Mat &A, Vec &u, Vec &rhs, Vec &u_exact, Vec &err,
        StructArgs &args, PetscLogEvent &warmUpEvent, PetscLogEvent &solvingEvent);

PetscErrorCode applyNeumannBC(Mat &A, Vec &RHS, const Vec &exact);
