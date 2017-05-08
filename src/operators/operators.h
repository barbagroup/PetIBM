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
# include "BodyPack.h"
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
PetscErrorCode createIdentity(const CartesianMesh &mesh, Mat &I);


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
 * \param bc an instance of Boundary class.
 * \param D returned divergence operator.
 * \param DCorrection returned matrix-free operator that serves as constraint mat.
 * \param normalize a bool indicating whether or not we want normalization 
 *        (default is True).
 *
 * Petsc matrix D should not be created before calling this function.
 *
 * \return PetscErrorCode.
 */
PetscErrorCode createDivergence(const CartesianMesh &mesh, const Boundary &bc, 
        Mat &D, Mat &DCorrection, 
        const PetscBool &normalize=PETSC_TRUE);


/**
 * \brief create an implicit Laplacian operator matrix.
 *
 * \param mesh an instance of CartesianMesh class.
 * \param bc an instance of Boundary class.
 * \param L returned matrix.
 * \param LCorrection returned matrix-free operator that serves as constraint mat.
 *
 * Petsc matrix L should not be created before calling this function.
 *
 * \return PetscErrorCode.
 */
PetscErrorCode createLaplacian(
        const CartesianMesh &mesh, const Boundary &bc, Mat &L, Mat &LCorrection);


/**
 * \brief create a matrix-free Mat for convective operater.
 *
 * \param mesh an instance of CartesianMesh.
 * \param bd an instance of Boundary.
 * \param H the returned matrix (an matrix-free Mat).
 *
 * Though this matrix-free operator is a PETSc Mat, not all PETSc Mat functions 
 * can be applied to it because we only register MatMult and MatDestroy 
 * functions in PetIBM. 
 *
 * \return PetscErrorCode.
 */
PetscErrorCode createConvection(
        const CartesianMesh &mesh, const Boundary &bd, Mat &H);


/**
 * \brief create non-normalized Bn.
 *
 * A is defined as \frac{I}{dt} - coeff \times Op, where I is identity matrix,
 * Op is the summed coefficient matrix for implicit operators, and coeff is the
 * coefficient from temporal integration. For example, for the case of Adams-
 * Bashforth for convective terms and Crank-Nicolson for diffusion term, 
 * AHead = \frac{I}{dt} - 0.5 \times LHead. If the implicit part is more complex,
 * there may be not only one coeff for each implicit operator, then, in order
 * to use this function, we may need to incorporate coefficients into each
 * implicit operator before calling this function. And then set coeff as 1.
 *
 * \param Op implicit operator.
 * \param dt size of a time step.
 * \param coeff coefficient of temporal integration for implicit operator.
 * \param n the order of Bn.
 * \param BnHead returned matrix.
 *
 * \return PetscErrorCode.
 */
PetscErrorCode createBnHead(const Mat &Op, const PetscReal &dt, 
        const PetscReal &coeff, const PetscInt &n, Mat &BnHead);


/**
 * \brief create normalized Bn.
 *
 * A is defined as \frac{inv(M)}{dt} - coeff \times Op, where M is defined as
 * MHead * inv(R), Op is the summed coefficient matrix for implicit operators, 
 * and coeff is the coefficient from temporal integration. For example, for the
 * case of Adams-Bashforth for convective terms and Crank-Nicolson for diffusion
 * term, A = \frac{inv(M)}{dt} - 0.5 \times L, where L = MHead * LHead * inv(R). 
 * If the implicit part is more complex,there may be not only one coeff for summed
 * implicit operator, then, in order to use this function, we may need to 
 * incorporate coefficients into each implicit operator before calling this 
 * function. And then set coeff as 1.
 *
 * \param Op implicit operator.
 * \param R the normalization matrix to define M = MHead * inv(R).
 * \param MHead the normalization matrix to define M = MHead * inv(R).
 * \param dt size of a time step.
 * \param coeff coefficient of temporal integration for implicit operator.
 * \param n the order of Bn.
 * \param Bn returned matrix.
 *
 * \return PetscErrorCode.
 */
PetscErrorCode createBn(const Mat &Op, const Mat &R, const Mat &MHead,
        const PetscReal &dt, const PetscReal &coeff, const PetscInt &n, 
        Mat &Bn);


/**
 * \brief create normalized Bn.
 *
 * A is defined as \frac{inv(M)}{dt} - coeff \times Op, where M is defined as
 * MHead * inv(R), Op is the summed coefficient matrix for implicit operators, 
 * and coeff is the coefficient from temporal integration. For example, for the
 * case of Adams-Bashforth for convective terms and Crank-Nicolson for diffusion
 * term, A = \frac{inv(M)}{dt} - 0.5 \times L, where L = MHead * LHead * inv(R). 
 * If the implicit part is more complex,there may be not only one coeff for summed
 * implicit operator, then, in order to use this function, we may need to 
 * incorporate coefficients into each implicit operator before calling this 
 * function. And then set coeff as 1.
 *
 * \param Op implicit operator.
 * \param M the first term in A = dt * inv(M) - coeff * Op.
 * \param dt size of a time step.
 * \param coeff coefficient of temporal integration for implicit operator.
 * \param n the order of Bn.
 * \param Bn returned matrix.
 *
 * \return PetscErrorCode.
 */
PetscErrorCode createBn(const Mat &Op, const Mat &M, const PetscReal &dt, 
        const PetscReal &coeff, const PetscInt &n, Mat &Bn);

/**
 * \brief create a Delta matrix.
 *
 * Rows represent Lagrangian points, while columes represent velocity mesh 
 * points.
 *
 * \param mesh an instance of CartesianMesh, which is background mesh.
 * \param bodies an instance of BodyPack.
 * \param Delta returned delta matrix.
 *
 * \return PetscErrorCode.
 */
PetscErrorCode createDelta(
        const CartesianMesh &mesh, const BodyPack &bodies, Mat &Delta);
