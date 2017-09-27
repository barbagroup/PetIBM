/*
 * singleboundaryneumann.h
 * Copyright (C) 2017 Pi-Yueh Chuang <pychuang@gwu.edu>
 *
 * Distributed under terms of the MIT license.
 */


# pragma once

// here goes headers from our PetIBM
# include <petibm/singleboundary.h>


namespace petibm
{
namespace boundary
{

class SingleBoundaryNeumann : public SingleBoundaryBase
{
public:

    SingleBoundaryNeumann(const type::Mesh &mesh, const type::BCLoc &loc, 
            const type::Field &field, const PetscReal &value); 

    virtual ~SingleBoundaryNeumann() = default;

protected:

    virtual PetscErrorCode setGhostICsKernel(
            const PetscReal &targetValue, type::GhostPointInfo &p);

    virtual PetscErrorCode updateEqsKernel(const PetscReal &targetValue,
            const PetscReal &dt, type::GhostPointInfo &p);

};

} // end of namespace boundary
} // end of namespace petibm
