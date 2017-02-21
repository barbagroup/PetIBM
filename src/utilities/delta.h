/*! Definition of the functions related to the delta functions.
 * \file delta.h
 */


#if !defined(DELTA_H)
#define DELTA_H

#include "types.h"

#include <petscsys.h>


// Returns the value of the discrete Delta function from Roma et al. (1999).
PetscReal dhRoma(PetscReal x, PetscReal h);

// Two-dimensional discrete delta function from Roma et al. (1999).
PetscReal delta(PetscReal x, PetscReal y, PetscReal hx, PetscReal hy);

// Three-dimensional discrete delta function from Roma et al. (1999).
PetscReal delta(PetscReal x, PetscReal y, PetscReal z, PetscReal hx, PetscReal hy, PetscReal hz);

// Defines if a point is in the disk or sphere of influence of another
// and calculates the displacement vector.
template<PetscInt dim>
PetscBool isInfluenced(PetscReal (&target)[dim], PetscReal (&source)[dim],
                       PetscReal (&maxDisp)[dim],
                       PetscReal (&widths)[dim], BoundaryType (&bTypes)[dim],
                       PetscReal *disp);

#endif
