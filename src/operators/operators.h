/***************************************************************************//**
 * \file operators.h
 * \author Anush Krishnan (anus@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Declarations of functions creating operators.
 */


// here goes PETSc headers
# include <petscmat.h>

// here goes headers from our PetIBM
# include "CartesianMesh.h"


/**
 * \brief create normalization matrix R.
 *
 * \param mesh an instance of CartesianMesh.
 * \param R returned matrix R.
 *
 * \return PetscErrorCode.
 */
PetscErrorCode createR(const CartesianMesh &mesh, Mat &R);
