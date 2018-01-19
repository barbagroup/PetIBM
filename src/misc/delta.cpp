/**
 * \file delta.cpp
 * \brief Implementations of descritized delta functions.
 * \author Anush Krishnan (anush@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \copyright MIT.
 */

// STL
# include <cmath>

// PetIBM
#include <petibm/delta.h>


namespace petibm
{
namespace delta
{

// 1D delta::Roma_et_al
PetscReal Roma_et_al(const PetscReal &rx, const PetscReal &drx)
{
    PetscReal r = std::abs(rx) / drx;

    if (r > 1.5)
        return 0.0;

    if (r > 0.5 && r <= 1.5)
        return (5.0 - 3.0 * r - 
                std::sqrt(-3.0 * (1.0 - r) * (1.0 - r) + 1.0)) / (6.0 * drx);

    return (1.0 + std::sqrt(-3.0 * r * r + 1.0)) / (3.0 * drx);
} // Roma_et_al


// 2D delta::Roma_et_al
PetscReal Roma_et_al(
        const PetscReal &rx, const PetscReal &drx,
        const PetscReal &ry, const PetscReal &dry)
{
    return delta::Roma_et_al(rx, drx) * delta::Roma_et_al(ry, dry);
} // Roma_et_al


// 3D delta::Roma_et_al(const PetscReal &, const PetscReal &, 
PetscReal Roma_et_al(
        const PetscReal &rx, const PetscReal &drx,
        const PetscReal &ry, const PetscReal &dry,
        const PetscReal &rz, const PetscReal &drz)
{
    return delta::Roma_et_al(rx, drx) * delta::Roma_et_al(ry, dry) * 
        delta::Roma_et_al(rz, drz);
} // Roma_et_al

} // end of namespace delta
} // end of namespace petibm
