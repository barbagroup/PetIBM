#if !defined(TYPES_H)
#define TYPES_H

/**
* @enum  BCType
* @brief Specifies the type of boundary condition
*/
enum BCType
{
	DIRICHLET,  ///< Dirichlet boundary condition
	NEUMANN,    ///< Neumann boundary condition
	CONVECTIVE, ///< Convective boundary condition
	PERIODIC    ///< Periodic boundary condition
};

/**
* @enum  boundary
* @brief Specifies the boundary of concern
*/
enum Boundary
{
	XMINUS,
	XPLUS,
	YMINUS,
	YPLUS,
	ZMINUS,
	ZPLUS
};

/**
* @enum  timeScheme
* @brief Numerical scheme used to discretise the time derivative
*/
enum TimeSteppingScheme
{
	EULER_EXPLICIT,    ///< Explicit Euler method (first order)
	EULER_IMPLICIT,    ///< Implicit Euler method (first order)
	ADAMS_BASHFORTH_2, ///< Second-order Adams-Bashforth scheme
	RUNGE_KUTTA_3,     ///< Third-order low storage Runge-Kutta method
	CRANK_NICOLSON     ///< Crank-Nicolson scheme (second order)
};

/**
* @enum  ibmScheme
* @brief The immersed boundary method used to solve for the flow.
*/
enum SolverType
{
	NAVIER_STOKES,  ///< No immersed bodies. Perot (1993)
	SAIKI_BIRINGEN, ///< Saiki & Biringen (1996)
	FADLUN_ET_AL,   ///< Fadlun et al (2000)
	TAIRA_COLONIUS, ///< Taira & Colonius (2007)
	SLL0,           ///< SLL0
	SLL1,           ///< SLL1
	SLL2            ///< SLL2
};

/**
* @enum  preconditionerType
* @brief Specify the type of preconditioner
*/
enum PreconditionerType
{
	NONE,                ///< No preconditioner
	DIAGONAL,            ///< Diagonal preconditioner
	SMOOTHED_AGGREGATION ///< Smoothed aggregation preconditioner
};

#endif