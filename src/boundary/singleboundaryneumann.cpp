/*
 * singleboundarydirichlet.cpp
 * Copyright (C) 2017 Pi-Yueh Chuang <pychuang@gwu.edu>
 *
 * Distributed under terms of the MIT license.
 */

#include <petibm/singleboundaryneumann.h>


namespace petibm
{
namespace boundary
{
    
    
SingleBoundaryNeumann::SingleBoundaryNeumann(
        const type::Mesh &inMesh, const type::BCLoc &inLoc,
        const type::Field &inField, const PetscReal &inValue):
    SingleBoundaryBase(inMesh, inLoc, inField, type::NEUMANN, inValue) {}


PetscErrorCode SingleBoundaryNeumann::setGhostICsKernel(
        const PetscReal &targetValue, type::GhostPointInfo &p)
{
    PetscFunctionBeginUser;
    
    PetscErrorCode  ierr;
    
    p.a0 = 1.0;
    p.a1 = normal * p.dL * value;
    
    p.value = p.a0 * targetValue + p.a1;
    
    PetscFunctionReturn(0);
}


PetscErrorCode SingleBoundaryNeumann::updateEqsKernel(const PetscReal &targetValue,
            const PetscReal &dt, type::GhostPointInfo &p)
{
    PetscFunctionBeginUser;
    // for time-independent Neumann BC, the coefficient a0 & a1 won't change
    PetscFunctionReturn(0);
}


}
}
