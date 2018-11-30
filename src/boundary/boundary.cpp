/**
 * \file boundary.cpp
 * \brief Implementation of boundary::BoundaryBase, type::Boundary, and factory
 * function. \copyright Copyright (c) 2016-2018, Barba group. All rights
 * reserved. \license BSD 3-Clause License.
 */

#include <petibm/boundary.h>
#include <petibm/boundarysimple.h>

namespace petibm
{
namespace boundary
{
// default destructor
BoundaryBase::~BoundaryBase()
{
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    PetscBool finalized;

    ierr = PetscFinalized(&finalized); CHKERRV(ierr);
    if (finalized) return;

    comm = MPI_COMM_NULL;
}  // ~BoundaryBase

PetscErrorCode BoundaryBase::destroy()
{
    PetscFunctionBeginUser;

    std::vector<std::vector<type::SingleBoundary>>().swap(bds);

    comm = MPI_COMM_NULL;
    mpiSize = mpiRank = 0;
    mesh.reset();

    PetscFunctionReturn(0);
}  // destroy

PetscErrorCode createBoundary(const type::Mesh &mesh, const YAML::Node &node,
                              type::Boundary &boundary)
{
    PetscFunctionBeginUser;

    boundary = std::make_shared<BoundarySimple>(mesh, node);

    PetscFunctionReturn(0);
}  // createBoundary

}  // end of namespace boundary
}  // end of namespace petibm
