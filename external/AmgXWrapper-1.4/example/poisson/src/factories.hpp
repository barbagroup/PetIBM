/**
 * \file factories.hpp
 * \brief definition of some functions generating something.
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \date 2017-06-01
 */


# pragma once

// PETSc
# include <petscsys.h>
# include <petscdmda.h>
# include <petscmat.h>
# include <petscvec.h>


/**
 * \brief generate coordinates of a 3D gird (PETSc DMDA) and cell sizes.
 *
 * \param grid [in] a PETSc DMDA.
 * \param Nx [in] number of cells in x direction.
 * \param Ny [in] number of cells in y direction.
 * \param Nz [in] number of cells in z direction.
 * \param Lx [in] length of domain in x direction.
 * \param Ly [in] length of domain in y direction.
 * \param Lz [in] length of domain in z direction.
 * \param dx [out] size of cells in x direction.
 * \param dy [out] size of cells in y direction.
 * \param dz [out] size of cells in z direction.
 * \param x [out] a PETSc Vec for coordinates in x direction.
 * \param y [out] a PETSc Vec for coordinates in y direction.
 * \param z [out] a PETSc Vec for coordinates in z direction.
 *
 * \return PetscErrorCode.
 */
PetscErrorCode generateGrid(const DM &grid, 
        const PetscInt &Nx, const PetscInt &Ny, const PetscInt &Nz,
        const PetscScalar &Lx, const PetscScalar &Ly, const PetscScalar &Lz,
        PetscScalar &dx, PetscScalar &dy, PetscScalar &dz,
        Vec &x, Vec &y, Vec &z);


/**
 * \brief generate coordinates of a 2D gird (PETSc DMDA) and cell sizes.
 *
 * \param grid [in] a PETSc DMDA.
 * \param Nx [in] number of cells in x direction.
 * \param Ny [in] number of cells in y direction.
 * \param Lx [in] length of domain in x direction.
 * \param Ly [in] length of domain in y direction.
 * \param dx [out] size of cells in x direction.
 * \param dy [out] size of cells in y direction.
 * \param x [out] a PETSc Vec for coordinates in x direction.
 * \param y [out] a PETSc Vec for coordinates in y direction.
 *
 * \return PetscErrorCode.
 */
PetscErrorCode generateGrid(const DM &grid, 
        const PetscInt &Nx, const PetscInt &Ny,
        const PetscScalar &Lx, const PetscScalar &Ly,
        PetscScalar &dx, PetscScalar &dy,
        Vec &x, Vec &y);


/**
 * \brief create a PETSc Vec for right-hand-side vector in 3D cases.
 *
 * \param grid [in] a PETSc DMDA.
 * \param x [in] a PETSc Vec for x coordinates.
 * \param y [in] a PETSc Vec for y coordinates.
 * \param z [in] a PETSc Vec for z coordinates.
 * \param rhs [out] a PETSc Vec for right-hand-side vector.
 *
 * \return PetscErrorCode.
 */
PetscErrorCode generateRHS(const DM &grid, 
        const Vec &x, const Vec &y, const Vec &z, Vec &rhs);


/**
 * \brief create a PETSc Vec for right-hand-side vector in 2D cases.
 *
 * \param grid [in] a PETSc DMDA.
 * \param x [in] a PETSc Vec for x coordinates.
 * \param y [in] a PETSc Vec for y coordinates.
 * \param rhs [out] a PETSc Vec for right-hand-side vector.
 *
 * \return PetscErrorCode.
 */
PetscErrorCode generateRHS(const DM &grid, 
        const Vec &x, const Vec &y, Vec &rhs);


/**
 * \brief generate a PETSc Vec for exact solution vector for 3D cases.
 *
 * \param grid [in] a PETSc DMDA.
 * \param x [in] a PETSc Vec for x coordinates.
 * \param y [in] a PETSc Vec for y coordinates.
 * \param z [in] a PETSc Vec for z coordinates.
 * \param exact [out] a PETSc Vec for exact solution vector.
 *
 * \return PetscErrorCode.
 */
PetscErrorCode generateExt(const DM &grid, 
        const Vec &x, const Vec &y, const Vec &z, Vec &exact);


/**
 * \brief generate a PETSc Vec for exact solution vector for 2D cases.
 *
 * \param grid [in] a PETSc DMDA.
 * \param x [in] a PETSc Vec for x coordinates.
 * \param y [in] a PETSc Vec for y coordinates.
 * \param exact [out] a PETSc Vec for exact solution vector.
 *
 * \return PetscErrorCode.
 */
PetscErrorCode generateExt(const DM &grid, 
        const Vec &x, const Vec &y, Vec &exact);


/**
 * \brief create a PETSc Mat for coefficient matrix for 3D cases.
 *
 * \param grid [in] a PETSc DMDA.
 * \param dx [in] cell size in x direction.
 * \param dy [in] cell size in y direction.
 * \param dz [in] cell size in z direction.
 * \param A [out] the PETSc Mat for coefficient matrix.
 *
 * \return  PetscErrorCode.
 */
PetscErrorCode generateA(const DM &grid, const PetscScalar &dx,
        const PetscScalar &dy, const PetscScalar &dz, Mat &A);


/**
 * \brief create a PETSc Mat for coefficient matrix for 2D cases.
 *
 * \param grid [in] a PETSc DMDA.
 * \param dx [in] cell size in x direction.
 * \param dy [in] cell size in y direction.
 * \param A [out] the PETSc Mat for coefficient matrix.
 *
 * \return  PetscErrorCode.
 */
PetscErrorCode generateA(const DM &grid, 
        const PetscScalar &dx, const PetscScalar &dy, Mat &A);
