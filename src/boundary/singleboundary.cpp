/**
 * \file singleboundary.cpp
 * \brief Implementation of the function `createSingleBoundary`.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 */

// here goes headers from our PetIBM
#include <petibm/singleboundary.h>
#include <petibm/singleboundaryconvective.h>
#include <petibm/singleboundarydirichlet.h>
#include <petibm/singleboundaryneumann.h>
#include <petibm/singleboundaryperiodic.h>

namespace petibm
{
namespace boundary
{
PetscErrorCode createSingleBoundary(const type::Mesh &mesh,
                                    const type::BCLoc &loc,
                                    const type::Field &field,
                                    const PetscReal &value,
                                    const type::BCType &bcType,
                                    type::SingleBoundary &singleBd)
{
    PetscFunctionBeginUser;

    switch (bcType)
    {
        case type::BCType::NOBC:
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                    "NOBC does not make sense here!");
            break;
        case type::BCType::PERIODIC:
            singleBd = std::make_shared<SingleBoundaryPeriodic>(mesh, loc,
                                                                field, value);
            break;
        case type::BCType::DIRICHLET:
            singleBd = std::make_shared<SingleBoundaryDirichlet>(mesh, loc,
                                                                 field, value);
            break;
        case type::BCType::NEUMANN:
            singleBd = std::make_shared<SingleBoundaryNeumann>(mesh, loc, field,
                                                               value);
            break;
        case type::BCType::CONVECTIVE:
            singleBd = std::make_shared<SingleBoundaryConvective>(mesh, loc,
                                                                  field, value);
            break;
    }

    PetscFunctionReturn(0);
}  // createSingleBoundary

}  // end of namespace boundary
}  // end of namespace petibm
