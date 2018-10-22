/**
 * \file lininterp.cpp
 * \brief Implementations of the interpolation classes and factory function.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 */

#include <petibm/lininterp.h>

namespace petibm
{

namespace misc
{
// Factory function to create a linear interpolation object.
PetscErrorCode createLinInterp(const MPI_Comm &comm,
                               const type::RealVec1D &point,
                               const type::Mesh &mesh,
                               const type::Field &field,
                               type::LinInterp &interp)
{
    PetscFunctionBeginUser;

    switch (mesh->dim)
    {
        case 2:
            interp = std::make_shared<BiLinInterp>(comm, point, mesh, field);
            break;
        case 3:
            interp = std::make_shared<TriLinInterp>(comm, point, mesh, field);
            break;
        default:
            SETERRQ(comm, PETSC_ERR_ARG_UNKNOWN_TYPE,
                    "Unknown number of dimensions. Accepted values are: 2 and 3");
    }

    PetscFunctionReturn(0);
}  // createLinInterp

//***************************************************************************//
//*************************      LinInterpBase      *************************//
//***************************************************************************//

// Destructor
LinInterpBase::~LinInterpBase()
{
    PetscErrorCode ierr;
    PetscBool finalized;

    PetscFunctionBeginUser;

    ierr = PetscFinalized(&finalized); CHKERRV(ierr);
    if (finalized) return;

    ierr = destroy(); CHKERRV(ierr);
}  // LinInterpBase::~LinInterpBase

// Manually destroy the data of the object.
PetscErrorCode LinInterpBase::destroy()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    comm = MPI_COMM_NULL;
    commSize = commRank = 0;

    ierr = VecDestroy(&coeffs); CHKERRQ(ierr);
    ierr = VecDestroy(&base); CHKERRQ(ierr);
    ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
    ierr = MatDestroy(&Op); CHKERRQ(ierr);
    ierr = VecDestroy(&sub); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // LinInterpBase::destroy

// Get the gridline indices of the front-bottom-left neighbor.
PetscErrorCode LinInterpBase::getBLGridlineIndices(const type::Mesh &mesh,
                                                   const type::Field &field)
{
    PetscFunctionBeginUser;

    std::vector<PetscReal>::iterator low;
    for (PetscInt d = 0; d < mesh->dim; ++d)
    {
        type::RealVec1D line(mesh->coord[field][d],
                             mesh->coord[field][d] + mesh->n[field][d]);
        low = std::lower_bound(line.begin(), line.end(), target[d]);
        idxDirs[d] = low - line.begin() - 1;
    }

    PetscFunctionReturn(0);
}  // LinInterpBase::getBLGridlineIndices

// Get the coordinates of the neighbors.
PetscErrorCode LinInterpBase::getBoxCoords(const type::Mesh &mesh,
                                           const type::Field &field)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = getBLGridlineIndices(mesh, field); CHKERRQ(ierr);

    for (PetscInt d = 0; d < mesh->dim; ++d)
    {
        bl[d] = mesh->coord[field][d][idxDirs[d]];
        tr[d] = mesh->coord[field][d][idxDirs[d] + 1];
    }

    PetscFunctionReturn(0);
}  // LinInterpBase::getBoxCoords

// Set up the direct solver to find the interpolation coefficients.
PetscErrorCode LinInterpBase::setUpKSP()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = KSPCreate(comm, &ksp); CHKERRQ(ierr);
    ierr = KSPSetType(ksp, KSPPREONLY); CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp, Op, Op); CHKERRQ(ierr);
    PC pc;
    ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);
    ierr = PCSetType(pc, PCLU); CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // LinInterpBase::setUpKSP

// Interpolate the field solution.
PetscErrorCode LinInterpBase::getValue(const DM &da, const Vec &vec, PetscReal &val)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = setSubVec(da, vec); CHKERRQ(ierr);
    ierr = KSPSolve(ksp, sub, coeffs); CHKERRQ(ierr);
    ierr = VecDot(coeffs, base, &val); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // LinInterpBase::getValue

//***************************************************************************//
//*************************      TriLinInterp       *************************//
//***************************************************************************//

// Constructor. Initialize the tri-linear interpolation object.
TriLinInterp::TriLinInterp(const MPI_Comm &comm,
                           const type::RealVec1D &point,
                           const type::Mesh &mesh,
                           const type::Field &field)
{
    init(comm, point, mesh, field);
}  // TriLinInterp::TriLinInterp

// Initialize the tri-linear interpolation object.
PetscErrorCode TriLinInterp::init(const MPI_Comm &inComm,
                                  const type::RealVec1D &point,
                                  const type::Mesh &mesh,
                                  const type::Field &field)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    comm = inComm;
    ierr = MPI_Comm_size(comm, &commSize); CHKERRQ(ierr);
    ierr = MPI_Comm_rank(comm, &commRank); CHKERRQ(ierr);

    target = point;
    idxDirs = type::IntVec1D(3, 0);
    bl = type::RealVec1D(3, 0.0);
    tr = type::RealVec1D(3, 0.0);
    base_a = type::RealVec1D(8, 0.0);
    Op_a = type::RealVec1D(64, 0.0);

    ierr = VecCreateSeq(comm, 8, &sub); CHKERRQ(ierr);
    ierr = VecDuplicate(sub, &base); CHKERRQ(ierr);
    ierr = VecDuplicate(sub, &coeffs);

    ierr = getBoxCoords(mesh, field); CHKERRQ(ierr);

    ierr = setUpBase(); CHKERRQ(ierr);
    ierr = setUpOp(); CHKERRQ(ierr);
    ierr = setUpKSP(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // TriLinInterp::init

// Set up the tri-linear base.
PetscErrorCode TriLinInterp::setUpBase()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    PetscReal x = target[0], y = target[1], z = target[2];
    base_a[0] = 1.0; base_a[1] = x; base_a[2] = y; base_a[3] = z;
    base_a[4] = x*y; base_a[5] = x*z; base_a[6] = y*z; base_a[7] = z*y*z;
    ierr = VecPlaceArray(base, &base_a[0]); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // TriLinInterp::setUpBase

// Set up the operator for a tri-linear interpolation.
PetscErrorCode TriLinInterp::setUpOp()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    PetscReal x0 = bl[0], y0 = bl[1], z0 = bl[2],
              x1 = tr[0], y1 = tr[1], z1 = tr[2];
    type::RealVec1D data = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                            x0, x1, x0, x1, x0, x1, x0, x1,
                            y0, y0, y1, y1, y0, y0, y1, y1,
                            z0, z0, z0, z0, z1, z1, z1, z1,
                            x0*y0, x1*y0, x0*y1, x1*y1, x0*y0, x1*y0, x0*y1, x1*y1,
                            x0*z0, x1*z0, x0*z0, x1*z0, x0*z1, x1*z1, x0*z1, x1*z1,
                            y0*z0, y0*z0, y1*z0, y1*z0, y0*z1, y0*z1, y1*z1, y1*z1,
                            x0*y0*z0, x1*y0*z0, x0*y1*z0, x1*y1*z0, x0*y0*z1, x1*y0*z1, x0*y1*z1, x1*y1*z1};
    Op_a = data;
    ierr = MatCreateSeqDense(comm, 8, 8, &Op_a[0], &Op); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // TriLinInterp::setUpOp

// Get the neighbor values.
PetscErrorCode TriLinInterp::setSubVec(const DM &da, const Vec &vec)
{
    PetscErrorCode ierr;
    PetscReal ***arr;
    PetscInt num = 0;

    PetscFunctionBeginUser;

    ierr = DMDAVecGetArrayRead(da, vec, &arr); CHKERRQ(ierr);
    PetscInt io = idxDirs[0], jo = idxDirs[1], ko = idxDirs[2];
    for (PetscInt k = 0; k < 2; ++k)
    {
        for (PetscInt j = 0; j < 2; ++j)
        {
            for (PetscInt i = 0; i < 2; ++i)
            {
                ierr = VecSetValue(sub, num, arr[ko + k][jo + j][io + i],
                                   INSERT_VALUES); CHKERRQ(ierr);
                num++;
            }
        }
    }
    ierr = DMDAVecRestoreArrayRead(da, vec, &arr); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // TriLinInterp::setSubVec

//***************************************************************************//
//*************************       BiLinInterp       *************************//
//***************************************************************************//

// Constructor. Initialize the bi-linear interpolation object.
BiLinInterp::BiLinInterp(const MPI_Comm &comm,
                         const type::RealVec1D &point,
                         const type::Mesh &mesh,
                         const type::Field &field)
{
    init(comm, point, mesh, field);
}  // BiLinInterp::BiLinInterp

// Initialize the bi-linear interpolation object.
PetscErrorCode BiLinInterp::init(const MPI_Comm &inComm,
                                 const type::RealVec1D &point,
                                 const type::Mesh &mesh,
                                 const type::Field &field)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    comm = inComm;
    ierr = MPI_Comm_size(comm, &commSize); CHKERRQ(ierr);
    ierr = MPI_Comm_rank(comm, &commRank); CHKERRQ(ierr);

    target = point;
    idxDirs = type::IntVec1D(2, 0);
    bl = type::RealVec1D(2, 0.0);
    tr = type::RealVec1D(2, 0.0);
    base_a = type::RealVec1D(4, 0.0);
    Op_a = type::RealVec1D(16, 0.0);

    ierr = VecCreateSeq(comm, 4, &sub); CHKERRQ(ierr);
    ierr = VecDuplicate(sub, &base); CHKERRQ(ierr);
    ierr = VecDuplicate(sub, &coeffs);

    ierr = getBoxCoords(mesh, field); CHKERRQ(ierr);

    ierr = setUpBase(); CHKERRQ(ierr);
    ierr = setUpOp(); CHKERRQ(ierr);
    ierr = setUpKSP(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // BiLinInterp::init

// Set up the bi-linear base.
PetscErrorCode BiLinInterp::setUpBase()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    PetscReal x = target[0], y = target[1];
    base_a[0] = 1.0; base_a[1] = x; base_a[2] = y; base_a[3] = x*y;
    ierr = VecPlaceArray(base, &base_a[0]); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // BiLinInterp::setUpBase

// Set up the operator for a bi-linear interpolation.
PetscErrorCode BiLinInterp::setUpOp()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    PetscReal x0 = bl[0], y0 = bl[1],
              x1 = tr[0], y1 = tr[1];
    type::RealVec1D data = {1.0, 1.0, 1.0, 1.0,
                            x0, x1, x0, x1,
                            y0, y0, y1, y1,
                            x0*y0, x1*y0, x1*y0, x1*y1};
    Op_a = data;
    ierr = MatCreateSeqDense(comm, 4, 4, &Op_a[0], &Op); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // BiLinInterp::setUpOp

// Get the neighbor values.
PetscErrorCode BiLinInterp::setSubVec(const DM &da, const Vec &vec)
{
    PetscErrorCode ierr;
    PetscReal **arr;
    PetscInt num = 0;

    PetscFunctionBeginUser;

    ierr = DMDAVecGetArrayRead(da, vec, &arr); CHKERRQ(ierr);
    PetscInt io = idxDirs[0], jo = idxDirs[1];
    for (PetscInt j = 0; j < 2; ++j)
    {
        for (PetscInt i = 0; i < 2; ++i)
        {
            ierr = VecSetValue(sub, num, arr[jo + j][io + i],
                               INSERT_VALUES); CHKERRQ(ierr);
            num++;
        }
    }
    ierr = DMDAVecRestoreArrayRead(da, vec, &arr); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // BiLinInterp::setSubVec

}  // end of namespace misc

}  // end of namespace petibm
