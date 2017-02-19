/*
 * misc.h
 * Copyright (C) 2017 Pi-Yueh Chuang <pychuang@gwu.edu>
 *
 * Distributed under terms of the MIT license.
 */

# pragma once

// here goes C++ STL
# include <vector>
# include <algorithm>
# include <cmath>

// here goes PETSc
# include <petscsys.h>

// here goes our own headers
# include "types.h"


/** \brief a namespace holding miscellaneous functions. */
namespace misc
{
    /** \brief find the pressure cell ID in which a given point is located (only in one direction).
     *
     * \param p the coordinate of the given point in a specific direction.
     * \param x the deviding points (vortex) of the pressure cells in a specific direction.
     * \param loc returned ID of the cell
     *
     * \return PetscErrorCode 
     *
     * \todo discuss if the point is in the first or the last cell, will PetIBM work?
     */
    PetscErrorCode findCell1D(
            const PetscReal &p, const std::vector<PetscReal> &x, PetscInt &loc);

    /** \brief  calculate and returned cell sizes of stretched grid in one direction.
     *
     * \param bg an input; start of the stretched region.
     * \param ed an input; end of the stretched region.
     * \param n an input; number of total cells in the stretched region.
     * \param r an input; value of the stretched ratio.
     * \param dL an output; cell sizes.
     *
     * \return PetscErrorCode
     */
    PetscErrorCode stretchGrid(
            const PetscReal &bg, const PetscReal &ed, 
            const PetscInt &n, const PetscReal &r, types::RealVec1D &dL);
}
