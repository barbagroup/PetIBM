/***************************************************************************//**
 * \file delta.h
 * \author Anush Krishnan (anus@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Definition of functions regarding to descritized delta functions.
 */


# pragma once


// PETSc
#include <petscsys.h>


namespace petibm
{
/** \brief namespace for all kinds of descritized delta functions. */
namespace delta
{
    /*!
     * \brief one-dimensional discrete delta function from Roma et al. (1999).
     *
     * \param rx distance between target and source.
     * \param drx window size.
     *
     * \returns The value of the discrete delta function.
     */
    PetscReal Roma_et_al(
            const PetscReal &rx, const PetscReal &drx);


    /*!
     * \brief two-dimensional discrete delta function from Roma et al. (1999).
     *
     * \param rx distance of the 1st component between target and source.
     * \param drx window size of the 1st component.
     * \param ry distance of the 2nd component between target and source.
     * \param dry window size of the 2nd component.
     *
     * \returns The value of the discrete delta function.
     */
    PetscReal Roma_et_al(
            const PetscReal &rx, const PetscReal &drx,
            const PetscReal &ry, const PetscReal &dry);


    /*!
     * \brief three-dimensional discrete delta function from Roma et al. (1999).
     *
     * \param rx distance of the 1st component between target and source.
     * \param drx window size of the 1st component.
     * \param ry distance of the 2nd component between target and source.
     * \param dry window size of the 2nd component.
     * \param rz distance of the 3rd component between target and source.
     * \param drz window size of the 3rd component.
     *
     * \returns The value of the discrete delta function.
     */
    PetscReal Roma_et_al(
            const PetscReal &rx, const PetscReal &drx,
            const PetscReal &ry, const PetscReal &dry,
            const PetscReal &rz, const PetscReal &drz);

} // end of namespace delta
} // end of namespace petibm
