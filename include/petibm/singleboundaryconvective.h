/***************************************************************************//**
 * \file singleboundary.h
 * \author Anush Krishnan (anus@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Definition of the class `SingleBoundary`.
 */


# pragma once

// here goes headers from our PetIBM
# include <petibm/singleboundary.h>


namespace petibm
{
namespace boundary
{

class SingleBoundaryConvective : public SingleBoundaryBase
{
public:

    SingleBoundaryConvective(const type::Mesh &mesh, const type::BCLoc &loc, 
            const type::Field &field, const PetscReal &value); 

    virtual ~SingleBoundaryConvective() = default;

protected:

    virtual PetscErrorCode setGhostICsKernel(
            const PetscReal &targetValue, type::GhostPointInfo &p);

    virtual PetscErrorCode updateEqsKernel(const PetscReal &targetValue,
            const PetscReal &dt, type::GhostPointInfo &p);

    PetscErrorCode (*kernel)(
            const PetscReal &targetValue, const PetscReal &dt, 
            const PetscReal &normal, const PetscReal &value,
            type::GhostPointInfo &p);
};

} // end of namespace boundary
} // end of namespace petibm
