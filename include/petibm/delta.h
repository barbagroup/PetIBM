/**
 * \file delta.h
 * \brief Prototype of Delta functions.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 */

#pragma once

#include <functional>
#include <vector>

#include <petscsys.h>

namespace petibm
{
/** \brief A namespace of all kinds of discretized delta functions.
 *
 * \ingroup miscModule
 */
namespace delta
{
/** \brief Regularized delta function from Roma et al. (1999).
 *
 * \param r [in] Distance between target and source
 * \param dr [in] Window size
 *
 * \returns The value of the regularized delta function.
 *
 * \ingroup miscModule
 */
PetscReal Roma_et_al_1999(const PetscReal &r, const PetscReal &dr);

/** \brief Regularized delta function from Peskin (2002).
 *
 * \param r [in] Distance between target and source
 * \param dr [in] Window size
 *
 * \returns The value of the regularized delta function.
 *
 * \ingroup miscModule
 */
PetscReal Peskin_2002(const PetscReal &r, const PetscReal &dr);

/** \brief Typedef to choose the regularized delta kernel to use. */
typedef std::function<PetscReal(const PetscReal &r,
                                const PetscReal &dr)> DeltaKernel;

/** \brief Get the delta kernel and size providing the name.
 *
 * \param name [in] Name of the delta kernel
 * \param kernel [out] The regularized delta kernel
 * \param size [out] Size of the delta kernel
*/
PetscErrorCode getKernel(const std::string &name,
                         DeltaKernel &kernel, PetscInt &size);

/** \brief Discrete delta function.
 *
 * \param source [in] Coordinates of the source point
 * \param target [in] Coordinates of the target point
 * \param h [in] Cell width
 *
 * \returns The value of the discrete delta function.
 *
 * \ingroup miscModule
 */
PetscReal delta(const std::vector<PetscReal> &source,
                const std::vector<PetscReal> &target,
                const std::vector<PetscReal> &widths,
                const DeltaKernel &kernel);

}  // end of namespace delta

}  // end of namespace petibm
