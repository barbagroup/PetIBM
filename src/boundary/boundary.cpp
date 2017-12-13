/*
 * boundary.cpp
 * Copyright (C) 2017 Pi-Yueh Chuang <pychuang@gwu.edu>
 *
 * Distributed under terms of the MIT license.
 */

# include <petibm/boundary.h>
# include <petibm/boundarysimple.h>


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
}


PetscErrorCode BoundaryBase::destroy()
{
    PetscFunctionBeginUser;
    PetscErrorCode ierr;

    std::vector<std::vector<type::SingleBoundary>>().swap(bds);

    comm = MPI_COMM_NULL;
    mpiSize = mpiRank = 0;
    mesh.reset();

    PetscFunctionReturn(0);
}

PetscErrorCode createBoundary(
        const type::Mesh &mesh, const YAML::Node &node,
        type::Boundary &boundary)
{
    PetscFunctionBeginUser;
    
    boundary = std::make_shared<BoundarySimple>(mesh, node);
    
    PetscFunctionReturn(0);
}

}
}
