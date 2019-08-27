/**
 * \file delta.cpp
 * \brief Implementations of descritized delta functions.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 */

#include <cmath>

#include <petibm/delta.h>

namespace petibm
{
namespace delta
{
// Regularized delta function from Roma et al. (1999).
PetscReal Roma_et_al_1999(const PetscReal &r, const PetscReal &dr)
{
    PetscReal x = std::abs(r) / dr;

    if (x > 1.5) return 0.0;

    if (x > 0.5 && x <= 1.5)
        return (5 - 3 * x - std::sqrt(-3 * (1 - x) * (1 - x) + 1)) / (6 * dr);

    return (1 + std::sqrt(-3 * x * x + 1)) / (3 * dr);
}  // Roma_et_al_1999

// Regularized delta function Peskin (2002).
PetscReal Peskin_2002(const PetscReal &r, const PetscReal &dr)
{
    PetscReal x = std::abs(r) / dr;

    if (x >= 0.0 && x <= 1.0)
        return (3 - 2 * x + std::sqrt(1 + 4 * x - 4 * x * x)) / (8 * dr);
    else if (x >= 1.0 && x <= 2.0)
        return (5 - 2 * x - std::sqrt(-7 + 12 * x - 4 * x * x)) / (8 * dr);
    return 0.0;
}  // Peskin_2002

// Get the delta kernel and size to use
PetscErrorCode getKernel(const std::string &name,
                         DeltaKernel &kernel, PetscInt &size)
{
    PetscFunctionBeginUser;

    if (name == "ROMA_ET_AL_1999")
    {
        kernel = Roma_et_al_1999;
        size = 2;
    }
    else if (name == "PESKIN_2002")
    {
        kernel = Peskin_2002;
        size = 3;
    }
    else
        SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_ARG_UNKNOWN_TYPE,
                 "No support for delta kernel `%s`.\n", name.c_str());

    PetscFunctionReturn(0);
}  // getKernel

// Discrete delta function.
PetscReal delta(const std::vector<PetscReal> &source,
                const std::vector<PetscReal> &target,
                const std::vector<PetscReal> &widths,
                const DeltaKernel &kernel)
{
    PetscReal phi = 1.0;
    for (unsigned int d = 0; d < widths.size(); ++d)
        phi *= kernel(source[d] - target[d], widths[d]);
    return phi;
}  // delta

}  // end of namespace delta
}  // end of namespace petibm
