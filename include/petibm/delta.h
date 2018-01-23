/**
 * \file delta.h
 * \brief Prototype of Delta functions.
 * \author Anush Krishnan (anush@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \copyright MIT.
 */


# pragma once


// PETSc
#include <petscsys.h>


namespace petibm
{
/**
 * \brief A namespace of all kinds of discretized delta functions.
 * \ingroup miscModule
 */
namespace delta
{
    /**
     * \brief One-dimensional discrete delta function from Roma et al. (1999).
     * \param rx [in] distance between target and source.
     * \param drx [in] window size.
     * \returns The value of the discrete delta function.
     * \ingroup miscModule
     */
    PetscReal Roma_et_al(
            const PetscReal &rx, const PetscReal &drx);


    /**
     * \brief Two-dimensional discrete delta function from Roma et al. (1999).
     * \param rx [in] distance of the 1st component between target and source.
     * \param drx [in] window size of the 1st component.
     * \param ry [in] distance of the 2nd component between target and source.
     * \param dry [in] window size of the 2nd component.
     * \returns The value of the discrete delta function.
     * \ingroup miscModule
     */
    PetscReal Roma_et_al(
            const PetscReal &rx, const PetscReal &drx,
            const PetscReal &ry, const PetscReal &dry);


    /**
     * \brief Three-dimensional discrete delta function from Roma et al. (1999).
     * \param rx [in] distance of the 1st component between target and source.
     * \param drx [in] window size of the 1st component.
     * \param ry [in] distance of the 2nd component between target and source.
     * \param dry [in] window size of the 2nd component.
     * \param rz [in] distance of the 3rd component between target and source.
     * \param drz [in] window size of the 3rd component.
     * \returns The value of the discrete delta function.
     * \ingroup miscModule
     */
    PetscReal Roma_et_al(
            const PetscReal &rx, const PetscReal &drx,
            const PetscReal &ry, const PetscReal &dry,
            const PetscReal &rz, const PetscReal &drz);

} // end of namespace delta
} // end of namespace petibm
