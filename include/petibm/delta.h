/**
 * \file delta.h
 * \brief Prototype of Delta functions.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 */

#pragma once

#include <petscsys.h>

namespace petibm
{
/**
 * \brief A namespace of all kinds of discretized delta functions.
 *
 * \ingroup miscModule
 */
namespace delta
{
/**
 * \brief One-dimensional discrete delta function from Roma et al. (1999).
 *
 * \param rx [in] Distance between target and source.
 * \param drx [in] Window size.
 *
 * \returns The value of the discrete delta function.
 *
 * \ingroup miscModule
 */
PetscReal Roma_et_al(const PetscReal &rx, const PetscReal &drx);

/**
 * \brief Two-dimensional discrete delta function from Roma et al. (1999).
 *
 * \param rx [in] Distance of the 1st component between target and source.
 * \param drx [in] Window size of the 1st component.
 * \param ry [in] Distance of the 2nd component between target and source.
 * \param dry [in] Window size of the 2nd component.
 *
 * \returns The value of the discrete delta function.
 *
 * \ingroup miscModule
 */
PetscReal Roma_et_al(const PetscReal &rx, const PetscReal &drx,
                     const PetscReal &ry, const PetscReal &dry);

/**
 * \brief Three-dimensional discrete delta function from Roma et al. (1999).
 *
 * \param rx [in] Distance of the 1st component between target and source.
 * \param drx [in] Window size of the 1st component.
 * \param ry [in] Distance of the 2nd component between target and source.
 * \param dry [in] Window size of the 2nd component.
 * \param rz [in] Distance of the 3rd component between target and source.
 * \param drz [in] Window size of the 3rd component.
 *
 * \returns The value of the discrete delta function.
 *
 * \ingroup miscModule
 */
PetscReal Roma_et_al(const PetscReal &rx, const PetscReal &drx,
                     const PetscReal &ry, const PetscReal &dry,
                     const PetscReal &rz, const PetscReal &drz);

}  // end of namespace delta

}  // end of namespace petibm
