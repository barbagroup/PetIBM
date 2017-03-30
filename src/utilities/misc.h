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
# include <functional>
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
    inline PetscErrorCode findCell1D(
            const PetscReal &p, const std::vector<PetscReal> &x, PetscInt &loc)
    {
        PetscFunctionBeginUser;

        if (x.back() < p || x[0] > p)
            SETERRQ1(PETSC_COMM_WORLD, 56, 
                    "body coordinate %e is outside domain !", p);

        auto it = std::lower_bound(x.begin(), x.end(), p);

        loc = it - x.begin() - 1;

        PetscFunctionReturn(0);
    }


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
    inline PetscErrorCode stretchGrid(
            const PetscReal &bg, const PetscReal &ed, 
            const PetscInt &n, const PetscReal &r, types::RealVec1D &dL)
    {
        PetscFunctionBeginUser;

        dL.resize(n);

        PetscReal       h;

        // calculate the size of the first cell
        dL[0] = (ed - bg) * (r - 1.0) / (std::pow(r, n) - 1.0);

        // dL[i] = dL[i-1] * r
        for(auto it=dL.begin()+1; it<dL.end(); ++it) *it = *(it -1) * r;

        PetscFunctionReturn(0);
    }


    /** \brief a helper struct to make looping function easier. */
    struct LoopBound {const PetscInt &bg, &ed;};


    /**
     * \brief a helper function to carry out a double loop on a given function.
     *
     * \param bound1 bound for the outer loop.
     * \param bound2 bound for the inner loop.
     * \param f the kernel that will be called in the function.
     *
     * \return PetscErrorCode.
     */
    inline PetscErrorCode doubleLoops(
            const LoopBound &bound1, const LoopBound &bound2, 
            const std::function<PetscErrorCode(const PetscInt &, const PetscInt &)> &f)
    {
        PetscFunctionBeginUser;

        PetscErrorCode      ierr;

        for(PetscInt i=bound1.bg; i<bound1.ed; i++)
            for(PetscInt j=bound2.bg; j<bound2.ed; j++)
            {
                ierr = f(i, j); CHKERRQ(ierr);
            }
                

        PetscFunctionReturn(0);
    }


    /**
     * \brief a helper function to carry out a triple loop on a given function.
     *
     * \param bound1 bound for the most outer loop.
     * \param bound2 bound for the middle loop.
     * \param bound3 bound for the most inner loop.
     * \param f the kernel that will be called in the function.
     *
     * \return PetscErrorCode.
     */
    inline PetscErrorCode tripleLoops(
            const LoopBound &bound1, const LoopBound &bound2,
            const LoopBound &bound3,
            const std::function<PetscErrorCode(
                const PetscInt &, const PetscInt &, const PetscInt &)> &f)
    {
        using namespace std::placeholders;

        PetscFunctionBeginUser;

        PetscErrorCode      ierr;

        for(PetscInt i=bound1.bg; i<bound1.ed; ++i)
        {
            ierr = doubleLoops(bound2, bound3, std::bind(f, i, _1, _2));
            CHKERRQ(ierr);
        }

        PetscFunctionReturn(0);
    }
}
