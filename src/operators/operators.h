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
# include "Boundary.h"
# include "types.h"


/**
 * \brief create normalization matrix R.
 *
 * \param mesh an instance of CartesianMesh.
 * \param R returned matrix R.
 *
 * Petsc matrix R should not be created before calling this function.
 *
 * \return PetscErrorCode.
 */
PetscErrorCode createR(const CartesianMesh &mesh, Mat &R);


/**
 * \brief create normalization matrix RInv.
 *
 * \param mesh an instance of CartesianMesh.
 * \param RInv returned matrix RInv.
 *
 * Petsc matrix RInv should not be created before calling this function.
 *
 * \return PetscErrorCode.
 */
PetscErrorCode createRInv(const CartesianMesh &mesh, Mat &RInv);


/**
 * \brief create normalization matrix MHead.
 *
 * \param mesh an instance of CartesianMesh.
 * \param MHead returned matrix MHead.
 *
 * Petsc matrix MHead should not be created before calling this function.
 *
 * \return PetscErrorCode.
 */
PetscErrorCode createMHead(const CartesianMesh &mesh, Mat &MHead);


/**
 * \brief create normalization matrix M.
 *
 * \param mesh an instance of CartesianMesh.
 * \param M returned matrix M.
 *
 * Petsc matrix M should not be created before calling this function.
 *
 * \return PetscErrorCode.
 */
PetscErrorCode createM(const CartesianMesh &mesh, Mat &M);


/**
 * \brief create identity matrix I.
 *
 * \param mesh an instance of CartesianMesh.
 * \param I returned matrix I.
 *
 * Petsc matrix I should not be created before calling this function.
 *
 * \return PetscErrorCode.
 */
PetscErrorCode createM(const CartesianMesh &mesh, Mat &I);


/**
 * \brief a function returning a gradient operator.
 *
 * \param mesh an instance of CartesianMesh.
 * \param G returned gradient operator.
 * \param normalize a bool indicating whether or not we want normalization 
 *        (default is True).
 *
 * Petsc matrix G should not be created before calling this function.
 *
 * \return PetscErrorCode.
 */
PetscErrorCode createGradient(const CartesianMesh &mesh, 
        Mat &G, const PetscBool &normalize=PETSC_TRUE);


/**
 * \brief a function returning a divergence operator.
 *
 * \param mesh an instance of CartesianMesh.
 * \param D returned divergence operator.
 * \param normalize a bool indicating whether or not we want normalization 
 *        (default is True).
 *
 * Petsc matrix D should not be created before calling this function.
 *
 * \return PetscErrorCode.
 */
PetscErrorCode createDivergence(const CartesianMesh &mesh, const Boundary &bc, 
        Mat &D, types::MatrixModifier &ghCoeff, 
        const PetscBool &normalize=PETSC_TRUE);


/**
 * \brief create an implicit Laplacian operator matrix.
 *
 * \param mesh an instance of CartesianMesh class.
 * \param bc an instance of Boundary class.
 * \param L returned matrix.
 *
 * Petsc matrix L should not be created before calling this function.
 *
 * \return PetscErrorCode.
 */
PetscErrorCode createLaplacian(
        const CartesianMesh &mesh, const Boundary &bc, 
        Mat &L, types::MatrixModifier &modifier);
