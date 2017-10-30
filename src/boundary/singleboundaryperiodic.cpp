/*
 * singleboundaryperiodic.cpp
 * Copyright (C) 2017 Pi-Yueh Chuang <pychuang@gwu.edu>
 *
 * Distributed under terms of the MIT license.
 */


#include <petibm/singleboundaryperiodic.h>


namespace petibm
{
namespace boundary
{
    
    
SingleBoundaryPeriodic::SingleBoundaryPeriodic(
        const type::Mesh &inMesh, const type::BCLoc &inLoc,
        const type::Field &inField, const PetscReal &inValue)
{
    // for Periodic, though we don't have to implement anything,
    // we still setup some basic values

    // create a shared pointer to mesh; bad practice...
    mesh = inMesh; 

    // obtain MPI information from CartesianMesh object
    comm = mesh->comm;
    mpiSize = mesh->mpiSize;
    mpiRank = mesh->mpiRank;

    // set dim
    dim = mesh->dim;

    // set the location
    loc = inLoc;
    
    // set field
    field = inField;
    
    // set type
    type = type::PERIODIC;
    
    // set value
    value = inValue;
    
    // set normal
    normal = ((int(loc)%2) == 0) ? -1.0 : 1.0;
}


PetscErrorCode SingleBoundaryPeriodic::setGhostICsKernel(
        const PetscReal &targetValue, type::GhostPointInfo &p)
{
    PetscFunctionBeginUser;
    // for Periodic BC, do nothing. Let PETSc handle that. 
    PetscFunctionReturn(0);
}


PetscErrorCode SingleBoundaryPeriodic::updateEqsKernel(const PetscReal &targetValue,
            const PetscReal &dt, type::GhostPointInfo &p)
{
    PetscFunctionBeginUser;
    // for Periodic BC, do nothing. Let PETSc handle that. 
    PetscFunctionReturn(0);
}


}
}
