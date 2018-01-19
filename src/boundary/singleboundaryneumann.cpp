/**
 * \file singleboundaryneumann.cpp
 * \brief Implementation of the class `SingleBoundaryNeumann`.
 * \author Anush Krishnan (anush@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \copyright MIT.
 */

#include <petibm/singleboundaryneumann.h>


namespace petibm
{
namespace boundary
{
    
    
SingleBoundaryNeumann::SingleBoundaryNeumann(
        const type::Mesh &inMesh, const type::BCLoc &inLoc,
        const type::Field &inField, const PetscReal &inValue):
    SingleBoundaryBase(inMesh, inLoc, inField, type::NEUMANN, inValue)
{

} // SingleBoundaryNeumann


PetscErrorCode SingleBoundaryNeumann::setGhostICsKernel(
        const PetscReal &targetValue, type::GhostPointInfo &p)
{
    PetscFunctionBeginUser;
    
    p.a0 = 1.0;
    p.a1 = normal * p.dL * value;
    
    p.value = p.a0 * targetValue + p.a1;
    
    PetscFunctionReturn(0);
} // setGhostICsKernel


PetscErrorCode SingleBoundaryNeumann::updateEqsKernel(const PetscReal &targetValue,
            const PetscReal &dt, type::GhostPointInfo &p)
{
    PetscFunctionBeginUser;
    // for time-independent Neumann BC, the coefficient a0 & a1 won't change
    PetscFunctionReturn(0);
} // updateEqsKernel

} // end of namespace boundary
} // end of namespace petibm
