/*
 * singleboundarydirichlet.cpp
 * Copyright (C) 2017 Pi-Yueh Chuang <pychuang@gwu.edu>
 *
 * Distributed under terms of the MIT license.
 */

#include <petibm/singleboundaryconvective.h>


PetscErrorCode kernelConvectiveDiffDir(
        const PetscReal &targetValue, const PetscReal &dt, 
        const PetscReal &normal, const PetscReal &value,
        petibm::type::GhostPointInfo &p)
{
    PetscFunctionBeginUser;
    
    p.a0 = -1.0;
    p.a1 = p.value + targetValue - 
        2.0 * normal * dt * value * (p.value - targetValue) / p.dL;
    
    PetscFunctionReturn(0);
}


PetscErrorCode kernelConvectiveSameDir(
        const PetscReal &targetValue, const PetscReal &dt, 
        const PetscReal &normal, const PetscReal &value,
        petibm::type::GhostPointInfo &p)
{
    PetscFunctionBeginUser;
    
    p.a0 = 0.0;
    p.a1 = p.value - normal * dt * value * (p.value - targetValue) / p.dL;
    
    PetscFunctionReturn(0);
}


namespace petibm
{
namespace boundary
{
    
    
SingleBoundaryConvective::SingleBoundaryConvective(
        const type::Mesh &inMesh, const type::BCLoc &inLoc,
        const type::Field &inField, const PetscReal &inValue):
    SingleBoundaryBase(inMesh, inLoc, inField, inValue)
{
    PetscInt dir = int(loc) / 2;
    
    if (dir == int(field))
        kernel = &kernelConvectiveSameDir;
    else
        kernel = &kernelConvectiveDiffDir;
}


PetscErrorCode SingleBoundaryConvective::setGhostICsKernel(
        const PetscReal &targetValue, type::GhostPointInfo &p)
{
    PetscFunctionBeginUser;
    
    PetscErrorCode  ierr;
    
    p.a1 = p.value = targetValue;
    
    PetscFunctionReturn(0);
}


PetscErrorCode SingleBoundaryConvective::updateEqsKernel(const PetscReal &targetValue,
            const PetscReal &dt, type::GhostPointInfo &p)
{
    PetscFunctionBeginUser;

    PetscErrorCode  ierr;
    
    ierr = kernel(targetValue, dt, normal, value, p); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}


}
}
