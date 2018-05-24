/*
 * \file private.h
 * \brief Prototypes of private functions that shared by different source files.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 */

# pragma once

# include <petscmat.h>

/**
 * \brief Operator< of MatStencil for using it as a key to map.
 *
 * Currently, this is made to be a private only used by PetIBM internally.
 *
 * For the detail of implementation, please refer to
 * std::lexicographical_compare. This is just a simplified version.
 */
static bool operator<(const MatStencil &l, const MatStencil &r)
{
    if (l.k < r.k) return true;
    if (l.k > r.k) return false;

    if (l.j < r.j) return true;
    if (l.j > r.j) return false;

    if (l.i < r.i) return true;
    if (l.i > r.i) return false;

    if (l.c < r.c) return true;

    return false;
} // operator<
