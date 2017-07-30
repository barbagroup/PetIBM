/***************************************************************************//**
 * \file delta.cpp
 * \author Anush Krishnan (anus@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief implementations of  descritized delta functions.
 */


// PetIBM
#include "delta.h"


namespace petibm
{
namespace utilities
{
namespace delta
{

/** \copydoc delta::Roma_et_al(const PetscReal &, const PetscReal &). */
PetscReal Roma_et_al(const PetscReal &rx, const PetscReal &drx)
{
    PetscReal r = fabs(rx) / drx;

    if (r > 1.5)
        return 0.0;

    if (r > 0.5 && r <= 1.5)
        return (5.0 - 3.0 * r - 
                sqrt(-3.0 * (1.0 - r) * (1.0 - r) + 1.0)) / (6.0 * drx);

    return (1.0 + sqrt(-3.0 * r * r + 1.0)) / (3.0 * drx);
} // Roma_et_al


/** \copydoc delta::Roma_et_al(const PetscReal &, const PetscReal &, 
 *           const PetscReal &, const PetscReal &). */
PetscReal Roma_et_al(
        const PetscReal &rx, const PetscReal &drx,
        const PetscReal &ry, const PetscReal &dry)
{
    return delta::Roma_et_al(rx, drx) * delta::Roma_et_al(ry, dry);
} // Roma_et_al


/** \copydoc delta::Roma_et_al(const PetscReal &, const PetscReal &, 
 *           const PetscReal &, const PetscReal &, const PetscReal &, 
 *           const PetscReal &). */
PetscReal Roma_et_al(
        const PetscReal &rx, const PetscReal &drx,
        const PetscReal &ry, const PetscReal &dry,
        const PetscReal &rz, const PetscReal &drz)
{
    return delta::Roma_et_al(rx, drx) * delta::Roma_et_al(ry, dry) * 
        delta::Roma_et_al(rz, drz);
} // Roma_et_al

} // end of namespace delta
} // end of namespace utilities
} // end of namespace petibm
