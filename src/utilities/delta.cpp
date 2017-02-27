/*! Implementation of the delta functions.
 * \file delta.cpp
 */

#include "delta.h"


/*!
 * \brief Returns the value of the discrete delta function
          from Roma et al. (1999).
 *
 * \param x Float at which the delta function is evaluated
 * \param h grid-spacing of the underlying Eulerian mesh
 *
 * \returns The value of the discrete delta function
 */
PetscReal dhRoma(PetscReal x, PetscReal h)
{
  PetscReal r = fabs(x)/h;
  if (r > 1.5)
    return 0.0;
  if (r > 0.5 && r <= 1.5)
    return 1.0/(6*h)*( 5.0 - 3.0*r - sqrt(-3.0*(1-r)*(1-r) + 1.0) );
  return 1.0/(3*h)*( 1.0 + sqrt(-3.0*r*r + 1.0) );
} // dhRoma


/*!
 * \brief Two-dimensional discrete delta function from Roma et al. (1999).
 *
 * \param x x-coordinate of the point at which the delta function is evaluated
 * \param y y-coordinate of the point at which the delta function is evaluated
 * \param hx grid-spacing of the underlying Eulerian mesh in the x-direction
 * \param hy grid-spacing of the underlying Eulerian mesh in the y-direction
 *
 * \returns The value of the 2D discrete delta function
 */
PetscReal delta(PetscReal x, PetscReal y, PetscReal hx, PetscReal hy)
{
  return dhRoma(x, hx) * dhRoma(y, hy);
} // delta


/*!
 * \brief Three-dimensional discrete delta function from Roma et al. (1999).
 *
 * \param x x-coordinate of the point at which the delta function is evaluated
 * \param y y-coordinate of the point at which the delta function is evaluated
 * \param z z-coordinate of the point at which the delta function is evaluated
 * \param hx grid-spacing of the underlying Eulerian mesh in the x-direction
 * \param hy grid-spacing of the underlying Eulerian mesh in the y-direction
 * \param hz grid-spacing of the underlying Eulerian mesh in the z-direction
 *
 * \returns The value of the 2D discrete delta function
 */
PetscReal delta(PetscReal x, PetscReal y, PetscReal z, PetscReal hx, PetscReal hy, PetscReal hz)
{
  return dhRoma(x, hx) * dhRoma(y, hy) * dhRoma(z, hz);
} // delta


/*! Defines if a point is in the disk or sphere of influence of another
    and calculates the displacement vector.
 *
 * \param target Coordinates of the target point
 * \param source Coordinates of the source point
 * \param maxDisp Dimensions of the domain of influence
 * \param widths Dimensions of the domain
 * \param bType Types of boundary conditions
 * \param disp The vector displacement to fill
 */
template<PetscInt dim>
PetscBool isInfluenced(PetscReal (&target)[dim], PetscReal (&source)[dim],
                       PetscReal (&maxDisp)[dim],
                       PetscReal (&widths)[dim], BoundaryType (&bTypes)[dim],
                       PetscReal *disp)
{
  PetscBool influenced = PETSC_TRUE;
  PetscInt i = 0;
  while (i < dim and influenced == PETSC_TRUE)
  {
    disp[i] = fabs(target[i] - source[i]);
    if (bTypes[i] == PERIODIC && disp[i] > widths[i]-disp[i])
      disp[i] = widths[i] - disp[i];
    if (disp[i] >= maxDisp[i])
      influenced = PETSC_FALSE;
    i++;
  }

  return influenced;
} // isInfluenced

template PetscBool isInfluenced<2>(PetscReal (&target)[2], PetscReal (&source)[2],
                                   PetscReal (&maxDisp)[2],
                                   PetscReal (&widths)[2], BoundaryType (&bTypes)[2],
                                   PetscReal *disp);
template PetscBool isInfluenced<3>(PetscReal (&target)[3], PetscReal (&source)[3],
                                   PetscReal (&maxDisp)[3],
                                   PetscReal (&widths)[3], BoundaryType (&bTypes)[3],
                                   PetscReal *disp);
