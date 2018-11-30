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
    PetscFunctionBeginUser;

    comm = MPI_COMM_NULL;
    commSize = commRank = 0;

    PetscFunctionReturn(0);
}  // LinInterpBase::destroy

// Initialize the interpolation object.
PetscErrorCode LinInterpBase::init(const MPI_Comm &inComm,
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
    idxDirs = type::IntVec1D(mesh->dim, 0);
    bl = type::RealVec1D(mesh->dim, 0.0);
    tr = type::RealVec1D(mesh->dim, 0.0);

    ierr = getBoxCoords(mesh, field); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // LinInterpBase::init

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

// Perform tri-linear interpolation.
PetscErrorCode TriLinInterp::interpolate(const DM &da, const Vec &vec, PetscReal &v)
{
    PetscErrorCode ierr;
    PetscReal ***a;
    PetscReal c00, c01, c10, c11, c0, c1;
    PetscInt &io = idxDirs[0], &jo = idxDirs[1], &ko = idxDirs[2];
    PetscReal &x = target[0], &y = target[1], &z = target[2];
    PetscReal &x0 = bl[0], &y0 = bl[1], &z0 = bl[2],
              &x1 = tr[0], &y1 = tr[1], &z1 = tr[2];
    PetscReal xd, yd, zd;

    PetscFunctionBeginUser;

    ierr = DMDAVecGetArrayRead(da, vec, &a); CHKERRQ(ierr);
    xd = (x - x0) / (x1 - x0);
    yd = (y - y0) / (y1 - y0);
    zd = (z - z0) / (z1 - z0);
    c00 = (1 - xd) * a[ko][jo][io] + xd * a[ko][jo][io+1];
    c01 = (1 - xd) * a[ko][jo+1][io] + xd * a[ko][jo+1][io+1];
    c10 = (1 - xd) * a[ko+1][jo][io] + xd * a[ko+1][jo][io+1];
    c11 = (1 - xd) * a[ko+1][jo+1][io] + xd * a[ko+1][jo+1][io+1];
    c0 = (1 - yd) * c00 + yd * c01;
    c1 = (1 - yd) * c10 + yd * c11;
    v = (1 - zd) * c0 + zd * c1;
    ierr = DMDAVecRestoreArrayRead(da, vec, &a); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // TriLinInterp::interpolate

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

// Perform bi-linear interpolation.
PetscErrorCode BiLinInterp::interpolate(const DM &da, const Vec &vec, PetscReal &v)
{
    PetscErrorCode ierr;
    PetscReal **a;
    PetscReal c0, c1;
    PetscInt &io = idxDirs[0], &jo = idxDirs[1];
    PetscReal &x = target[0], &y = target[1];
    PetscReal &x0 = bl[0], &y0 = bl[1],
              &x1 = tr[0], &y1 = tr[1];
    PetscReal xd, yd;

    PetscFunctionBeginUser;

    ierr = DMDAVecGetArrayRead(da, vec, &a); CHKERRQ(ierr);
    xd = (x - x0) / (x1 - x0);
    yd = (y - y0) / (y1 - y0);
    c0 = (1 - xd) * a[jo][io] + xd * a[jo][io + 1];
    c1 = (1 - xd) * a[jo + 1][io] + xd * a[jo + 1][io + 1];
    v = (1 - yd) * c0 + yd * c1;
    ierr = DMDAVecRestoreArrayRead(da, vec, &a); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // BiLinInterp::interpolate

}  // end of namespace misc

}  // end of namespace petibm
