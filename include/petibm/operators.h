/**
 * \file operators.h
 * \brief Prototypes of factory functions for operators.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 */


// here goes PETSc headers
# include <petscmat.h>

// here goes headers from our PetIBM
# include <petibm/mesh.h>
# include <petibm/boundary.h>
# include <petibm/bodypack.h>
# include <petibm/type.h>


/**
 * \defgroup operatorModule Operator factories
 * \brief Matrices & Operators used in Perot's framework (1993, 2002).
 * 
 * An operator is a matrix, which will be the type of PETSc Mat. It may be a 
 * sparse matrix or a matrix-free matrix. Operators in CFD are based on the 
 * discretization of domain and boundary conditions, i.e. mesh. Therefore, all
 * operator factories require a \ref petibm::type::Mesh "Mesh" instance as an 
 * input, and most of them require a \ref petibm::type::Boundary "Boundary"
 * instance as another input.
 * 
 * Given that PetIBM now only supports staggered grid, most of the operators
 * only make sense to staggered Cartesian CFD methods. So the operator factories
 * can not create the operators for arbitrary fields. For example, 
 * \ref petibm::operators::createLaplacian "createLaplacian" can only create 
 * Laplacian for velocity fields because in our CFD framework (Perot 1993), only 
 * the Laplacian of velocity fields is required. This factory function cannot 
 * create Laplacian operators for arbitrary fields.
 * 
 * For operators taking boundary conditions, we decouple the interior part and
 * boundary corrections. For example, applying a divergence operator to velocity
 * will become \f[(D + D_{bc})u\f] in which \f$D\f$ and \f$D_{bc}\f$ represent
 * divergence operators of interior points and ghost points (i.e., boundary 
 * correction), respectively.
 * 
 * Boundary conditions in PetIBM are implemented with ghost point schemes. Most
 * of boundary conditions can be represented with, for example, 
 * \f[ u_{i} = a_{0} u_{i-1} + a_{1} \f]
 * where \f$u_{i}\f$ and \f$u_{i-1}\f$ represent a ghost point and its 
 * corresponding interior point, respectively. Hence when we apply ghost point 
 * schemes to operators, the \f$a_{0}\f$ part will get into the interior operator 
 * (e.g. \f$D\f$), and the \f$a_{1}\f$ part will be decoupled into the boundary 
 * correction matrix, i.e., \f$D_{bc}\f$.
 * 
 * Here is an example code for creating a divergence operator and its boundary 
 * correction operator, and calculating the divergence of velocity.
 * 
 * \code
 * PetscErrorCode ierr;
 * Mat D, Dbc;
 * Vec div;
 * petibm::type::Solution solution;
 * petibm::type::Mesh mesh;
 * petibm::type::Boundary boundary;
 * petibm::type::Solution solution;
 * 
 * // create mesh with petibm::mesh::createMesh
 * // create boundary with petibm::boundary::createBoundary
 * // create solution with petibm::solution::createSolution, and have some values 
 * // create vector div with correct memory allocation
 * 
 * // create divergence operator
 * ierr = petibm::createDivergence(mesh, boundary, D, Dbc); CHKERRQ(ierr);
 * 
 * // calculate divergence of velocity, result will be in div
 * ierr = MatMult(D, solution->UGlobal, div); CHKERRQ(ierr);
 * ierr = MatMultAdd(Dbc, solution->UGlobal, div, div); CHKERRQ(ierr);
 * \endcode
 * 
 * \ingroup petibm
 */


namespace petibm
{
/**
 * \brief Collections of factory functions of operators.
 * \ingroup operatorModule
 */
namespace operators
{

/**
 * \brief Create a matrix of flux areas, \f$R\f$.
 * \param mesh [in] an instance of petibm::type::Mesh.
 * \param R [out] a PETSc Mat; returned matrix \f$R\f$.
 * \return PetscErrorCode.
 *
 * PETSc matrix R should not be created before calling this function.
 * 
 * Operator \f$R\f$ represents the flux areas of velocity points. For
 * example, if \f$q\f$ is vector of velocity flux, and \f$u\f$ is velocity 
 * vector, then \f$Rq=u\f$.
 * 
 * \ingroup operatorModule
 * \see petibm::operators::createRInv
 */
PetscErrorCode createR(const type::Mesh &mesh, Mat &R);


/**
 * \brief Create a matrix of inversed flux areas, \f$R^{-1}\f$.
 * \param mesh [in] an instance of petibm::type::Mesh.
 * \param RInv [out] a PETSc Mat; returned matrix \f$R^{-1}\f$.
 * \return PetscErrorCode.
 *
 * PETSc matrix RInv should not be created before calling this function.
 * 
 * Operator \f$R^{-1}\f$ is the inverse of \f$R\f$. And it has relationship
 * \f$R^{-1}u=q\f$, in which \f$u\f$ and \f$q\f$ are as defined in 
 * petibm::operators::createR.
 * 
 * \ingroup operatorModule
 * \see petibm::operators::createR
 */
PetscErrorCode createRInv(const type::Mesh &mesh, Mat &RInv);


/**
 * \brief Create a matrix of cell widths of velocity points, \f$\hat{M}\f$.
 * \param mesh [in] an instance of petibm::type::Mesh.
 * \param MHead [out] a PETSc Mat; returned matrix \f$\hat{M}\f$.
 * \return PetscErrorCode.
 *
 * PETSc matrix MHead should not be created before calling this function.
 * 
 * Operator \f$\hat{M}\f$ represents the cell width of velocity points. The
 * cell width of a \f$u\f$-velocity point with index \f$(i, j, k)\f$, for 
 * example, is \f$dx_{i,j,k}\f$, and the cell width of \f$v\f$-velocity at
 * \f$(i, j, k)\f$ is \f$dy_{i,j,k}\f$.
 * 
 * The multiplication of \f$\hat{M}\f$ and \f$R\f$ returns the cell volumes of 
 * velocity points, i.e., \f$\hat{M}R=cell\ volumes\f$.
 *
 * \ingroup operatorModule
 */
PetscErrorCode createMHead(const type::Mesh &mesh, Mat &MHead);


/**
 * \brief Create a matrix \f$M=\hat{M}R^{-1}\f$.
 * \param mesh [in] an instance of petibm::type::Mesh.
 * \param M [out] a PETSc Mat; returned matrix \f$M\f$.
 * \return PetscErrorCode.
 *
 * PETSc matrix M should not be created before calling this function.
 * 
 * In some CFD schemes, the combination of \f$\hat{M}\f$ and \f$R^{-1}\f$
 * is used frequently. So this function creates a matrix \f$M=\hat{M}R^{-1}\f$
 * for convenience.
 *
 * \ingroup operatorModule
 */
PetscErrorCode createM(const type::Mesh &mesh, Mat &M);


/**
 * \brief Create an identity matrix \f$I\f$ for velocity fields.
 * \param mesh [in] an instance of petibm::type::Mesh.
 * \param I [out] a PETSc Mat; returned matrix \f$I\f$.
 * \return PetscErrorCode.
 *
 * PETSc matrix I should not be created before calling this function.
 *
 * \ingroup operatorModule
 */
PetscErrorCode createIdentity(const type::Mesh &mesh, Mat &I);


/**
 * \brief Create a gradient operator, \f$G\f$, for pressure field.
 * \param mesh [in] an instance of petibm::type::Mesh.
 * \param G [out] a PETSc Mat; returned gradient operator, \f$G\f$.
 * \param normalize [in] an optional PetscBool indicating whether or not we want 
 *                       normalization (default is PETSC_TRUE).
 * \return PetscErrorCode.
 *
 * Petsc matrix G should not be created before calling this function.
 * 
 * If `normalize=PETSC_TRUE` (which is default), then the gradient operator, 
 * \f$G\f$, will only contain `0` and `1`.
 * 
 * This gradient operator doesn't require boundary conditions because we only
 * need gradient of pressure at interior velocity points (see Perot 1993). And 
 * it can only be applied to pressure field. For example:
 * \code
 * ierr = petibm::operators::createGradient(mesh, G); CHKERRQ(ierr);
 * ierr = MatMult(G, solution->pGlobal, grad); CHKERRQ(ierr);
 * \endcode
 *
 * \ingroup operatorModule
 */
PetscErrorCode createGradient(const type::Mesh &mesh, 
                              Mat &G,
                              const PetscBool &normalize=PETSC_TRUE);


/**
 * \brief Create a divergence operator, \f$D\f$, and corresponding boundary
 *        correction, \f$D_{bc}\f$, for velocity fields.
 * \param mesh [in] an instance of petibm::type::Mesh.
 * \param bc [in] an instance of petibm::type::Boundary.
 * \param D [out] a PETSc Mat; returned divergence operator, \f$D\f$.
 * \param DCorrection [out] a PETSc Mat; returned boundary correction, \f$D_{bc}\f$.
 * \param normalize [in] an optional PetscBool indicating whether or not we want 
 *                  normalization (default is PETSC_TRUE).
 * \return PetscErrorCode.
 *
 * PETSc matrix D should not be created before calling this function.
 * 
 * If `normalize=PETSC_TRUE` (which is default), then the divergence operator, 
 * \f$D\f$, will only contain `0` and `1`.
 * 
 * The divergence operator calculates the divergence of velocity on pressure 
 * points. The complete divergence operator is \f$(D+D_{bc})\f$. However, the 
 * boundary correction, \f$D_{bc}\f$, represents the constant coefficients (it 
 * means the coefficients that have nothing to do with interior points, though
 * it may still be time-dependent) in ghost point schemes, so it's safe to put
 * the boundary correction to right-hand side in a linear system.
 *
 * \ingroup operatorModule
 */
PetscErrorCode createDivergence(const type::Mesh &mesh,
                                const type::Boundary &bc, 
                                Mat &D,
                                Mat &DCorrection, 
                                const PetscBool &normalize=PETSC_TRUE);


/**
 * \brief Create a Laplacian operator, \f$L\f$, and boundary correction, 
 *        \f$L_{bc}\f$ for velocity fields.
 * \param mesh [in] an instance of petibm::type::Mesh.
 * \param bc [in] an instance of petibm::type::Boundary.
 * \param L [out] a PETSc Mat; returned Laplacian, \f$L\f$.
 * \param LCorrection [out] a PETSc Mat; returned boundary correction, \f$L_{bc}\f$.
 * \return PetscErrorCode.
 *
 * PETSc matrix L should not be created before calling this function.
 * 
 * The complete Laplacian is \f$(L+L_{bc})\f$. The boundary correction represents 
 * the constant coefficients (the coefficients that have nothing to do with 
 * interior velocity points) in ghost point schemes. So it's often seen the 
 * boundary correction \f$L_{bc}\f$ on right-hand side of a system, while \f$L\f$
 * on left-hand side.
 *
 * \ingroup operatorModule
 */
PetscErrorCode createLaplacian(const type::Mesh &mesh,
                               const type::Boundary &bc,
                               Mat &L,
                               Mat &LCorrection);


/**
 * \brief Create a matrix-free Mat for convection operator, \f$H\f$.
 * \param mesh [in] an instance of petibm::type::Mesh.
 * \param bc [in] an instance of petibm::type::Boundary.
 * \param H [out] a PETSc Mat; the returned operator, \f$H\f$.
 * \return PetscErrorCode.
 *
 * Though this matrix-free operator is a PETSc Mat, not all PETSc Mat functions 
 * can be applied to it because we only register MatMult and MatDestroy 
 * functions in PetIBM. 
 * 
 * If \f$u\f$ is the composite velocity vector, then the multiplication \f$Hu\f$
 * will return the convection terms at velocity points.
 *
 * \ingroup operatorModule
 */
PetscErrorCode createConvection(const type::Mesh &mesh,
                                const type::Boundary &bc,
                                Mat &H);


/**
 * \brief Create non-normalized matrix of approximated \f$A^{-1}\f$, \f$B_n\f$.
 * \param Op [in] a PETSc Mat; implicit operator.
 * \param dt [in] size of a time step.
 * \param coeff [in] coefficient of temporal integration for implicit operator.
 * \param n [in] integer; the order of Bn.
 * \param BnHead [out] returned \f$B_n\f$.
 * \return PetscErrorCode.
 *
 * \f$A\f$ is defined as * \f[A=\frac{I}{dt} - coeff \times Op\f] 
 * where \f$I\f$ is an identity matrix, \f$Op\f$ is the summation of all 
 * implicit operators, and \f$coeff\f$ is the coefficient from temporal integration. 
 * For example, for the case of Adams-Bashforth for convective terms and 
 * Crank-Nicolson for diffusion term, 
 * \f[A = \frac{I}{dt} - 0.5 \times L\f]
 * 
 * If the implicit part is more complicated, there may be not only one coefficient 
 * for each implicit operator. Then, in order to use this function, we may need to 
 * incorporate coefficients and all implicit operators into a big implicit operator 
 * before calling this function. And then set \f$coeff\f$ as \f$1\f$. For example, 
 * if an \f$A\f$ is defined as
 * \f[ A = \frac{I}{dt} - c_{1}\times Op_{1} - c_{2}\times Op_{2}\f]
 * then we need to first define \f[Op = c_1\times Op_1 + c_2\times Op_2\f]
 * so that \f$A\f$ can have the form of
 * \f[ A = \frac{I}{dt} - 1.0\times Op \f].
 * 
 * \ingroup operatorModule
 */
PetscErrorCode createBnHead(const Mat &Op, const PetscReal &dt, 
        const PetscReal &coeff, const PetscInt &n, Mat &BnHead);


/**
 * \brief Create normalized matrix of approximated \f$A^{-1}\f$, \f$B_n\f$.
 * \param Op [in] a PETSc Mat; implicit operator.
 * \param R [in] the normalization matrix to define \f$M = \hat{M}R^{-1}\f$.
 * \param MHead [in] the normalization matrix to define \f$M = \hat{M}R^{-1}\f$.
 * \param dt [in] size of a time step.
 * \param coeff [in] coefficient of temporal integration for implicit operator.
 * \param n [in] an integer; the order of \f$B_n\f$.
 * \param Bn [out] a PETSc Mat; returned matrix, \f$B_n\f$.
 * \return PetscErrorCode.
 *
 * Similar to \ref petibm::operators::createBnHead "createBnHead", except that
 * A is now defined as 
 * \f[ A = \frac{M}{dt} - coeff \times \hat{M}OpR \f]
 * where 
 * \f[ M = \hat{M}R^{-1} \f] 
 * 
 * \ingroup operatorModule
 */
PetscErrorCode createBn(const Mat &Op, const Mat &R, const Mat &MHead,
        const PetscReal &dt, const PetscReal &coeff, const PetscInt &n, 
        Mat &Bn);


/**
 * \brief Create normalized matrix of approximated \f$A^{-1}\f$, \f$B_n\f$.
 * \param Op [in] a PETSc Mat; implicit operator.
 * \param M [in] the \f$M\f$ matrix as defined in \f$A = M/dt - coeff * Op\f$.
 * \param dt [in] size of a time step.
 * \param coeff [in] coefficient of temporal integration for implicit operator.
 * \param n [in] an integer; the order of Bn.
 * \param Bn [out] returned matrix, \f$B_n\f$.
 * \return PetscErrorCode.
 *
 * Similar to \ref petibm::operators::createBnHead "createBnHead", except that
 * A is now defined as \f[ A = \frac{M}{dt} - coeff \times Op \f]
 *
 * \ingroup operatorModule
 */
PetscErrorCode createBn(const Mat &Op, const Mat &M, const PetscReal &dt, 
        const PetscReal &coeff, const PetscInt &n, Mat &Bn);

/**
 * \brief Create a Delta operator, \f$Delta\f$.
 * \param mesh [in] an instance of petibm::type::Mesh.
 * \param bc [in] an instance of petibm::type::Boundary.
 * \param bodies [in] an instance of petibm::type::BodyPack.
 * \param Delta [out] a PETSc Mat; returned \f$Delta\f$.
 * \return PetscErrorCode.
 *
 * In \f$Delta\f$, rows represent Lagrangian points, while columns represent 
 * velocity mesh points.
 * 
 * An entry in \f$Delta\f$, says \f$Delta_{(i,j)}\f$, is defined as
 * \f[ Delta_{(i,j)} = d(x_j-\xi_i)d(y_j-\eta_i)d(z_j-\zeta_i) \f]
 * where \f$(x_j,y_j,z_j)\f$ and \f$(\xi_i,\eta_i,\zeta_i)\f$ represent the
 * coordinates of the \f$j\f$-th Eulerian grid point and the \f$i\f$-th 
 * Lagrangian point, respectively.
 * \f$d(...)\f$ is discretized delta function following Roma et. al. (1999).
 *
 * \ingroup operatorModule
 */
PetscErrorCode createDelta(const type::Mesh &mesh,
                           const type::Boundary &bc,
                           const type::BodyPack &bodies,
                           Mat &Delta);

} // end of namespace operators
} // end of namespace petibm
