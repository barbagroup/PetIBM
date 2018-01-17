/**
 * \file singleboundarydirichlet.cpp
 * \brief Implementation of the class `SingleBoundaryDirichlet`.
 * \author Anush Krishnan (anus@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \copyright MIT.
 */

#include <petibm/singleboundarydirichlet.h>


namespace petibm
{
namespace boundary
{
    
    
SingleBoundaryDirichlet::SingleBoundaryDirichlet(
        const type::Mesh &inMesh, const type::BCLoc &inLoc,
        const type::Field &inField, const PetscReal &inValue):
    SingleBoundaryBase(inMesh, inLoc, inField, type::DIRICHLET, inValue) {}


PetscErrorCode SingleBoundaryDirichlet::setGhostICsKernel(
        const PetscReal &targetValue, type::GhostPointInfo &p)
{
    PetscFunctionBeginUser;
    
    PetscErrorCode  ierr;
    
    PetscInt    dir = int(loc) / 2;
    
    // bad idea; this means every ghost point, p, passed in will need to go
    // through this if-condition, even though they all have same dir and field
    // Luckly, for time-independent Dirichlet BC, a0 & a1 won't change afterward
    if (dir == int(field))
    {    
        p.a0 = 0.0;
        p.a1 = value;
    }
    else
    {
        p.a0 = - 1.0;
        p.a1 = 2.0 * value;
    }
    
    p.value = p.a0 * targetValue + p.a1;
    
    PetscFunctionReturn(0);
}


PetscErrorCode SingleBoundaryDirichlet::updateEqsKernel(const PetscReal &targetValue,
            const PetscReal &dt, type::GhostPointInfo &p)
{
    PetscFunctionBeginUser;
    // for time-independent Dirichlet BC, the coefficient a0 & a1 won't change
    PetscFunctionReturn(0);
}


}
}
