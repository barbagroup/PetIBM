/**
 * \file singleboundaryconvective.cpp
 * \brief Implementation of the class `SingleBoundaryConvective`.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 */

#include <petibm/singleboundaryconvective.h>

namespace  // anonymous namespace for internal linkage
{
PetscErrorCode kernelConvectiveDiffDir(const PetscReal &targetValue,
                                       const PetscReal &dt,
                                       const PetscReal &normal,
                                       const PetscReal &value,
                                       petibm::type::GhostPointInfo &p)
{
    PetscFunctionBeginUser;

    p.a0 = -1.0;
    p.a1 = p.value + targetValue -
           2.0 * normal * dt * value * (p.value - targetValue) / p.dL;

    PetscFunctionReturn(0);
}  // kernelConvectiveDiffDir

PetscErrorCode kernelConvectiveSameDir(const PetscReal &targetValue,
                                       const PetscReal &dt,
                                       const PetscReal &normal,
                                       const PetscReal &value,
                                       petibm::type::GhostPointInfo &p)
{
    PetscFunctionBeginUser;

    p.a0 = 0.0;
    p.a1 = p.value - normal * dt * value * (p.value - targetValue) / p.dL;

    PetscFunctionReturn(0);
}  // kernelConvectiveSameDir
}  // end of anonymous namespace

namespace petibm
{
namespace boundary
{
SingleBoundaryConvective::SingleBoundaryConvective(const type::Mesh &inMesh,
                                                   const type::BCLoc &inLoc,
                                                   const type::Field &inField,
                                                   const PetscReal &inValue)
    : SingleBoundaryBase(inMesh, inLoc, inField, type::CONVECTIVE, inValue)
{
    PetscInt dir = int(loc) / 2;

    if (dir == int(field))
        kernel = &kernelConvectiveSameDir;
    else
        kernel = &kernelConvectiveDiffDir;
}  // SingleBoundaryConvective

PetscErrorCode SingleBoundaryConvective::destroy()
{
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    kernel = nullptr;
    ierr = SingleBoundaryBase::destroy(); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}  // destroy

PetscErrorCode SingleBoundaryConvective::setGhostICsKernel(
    const PetscReal &targetValue, type::GhostPointInfo &p)
{
    PetscFunctionBeginUser;

    PetscErrorCode ierr;

    // at beginning (t=0), due to lack of information of previous solution,
    // we just assume the ghost point has the same value as target boundary
    // point does
    p.value = targetValue;

    // due to p.value=targetValue, it doesn't matter how big is time-step size
    ierr = kernel(targetValue, 0.0, normal, value, p); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // setGhostICsKernel

PetscErrorCode SingleBoundaryConvective::updateEqsKernel(
    const PetscReal &targetValue, const PetscReal &dt, type::GhostPointInfo &p)
{
    PetscFunctionBeginUser;

    PetscErrorCode ierr;

    ierr = kernel(targetValue, dt, normal, value, p); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // updateEqsKernel

}  // end of namespace boundary
}  // end of namespace petibm
