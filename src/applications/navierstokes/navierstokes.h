/**
 * \file navierstokes.h
 * \brief Definition of the class \c NavierStokesSolver.
 */

#pragma once

#include "utilities/CartesianMesh.h"
#include "utilities/FlowDescription.h"
#include "utilities/SimulationParameters.h"
#include "utilities/Solutions.h"
#include "utilities/Boundary.h"
#include "utilities/TimeIntegration.h"
#include "linSolvers/LinSolver.h"
#include "operators/operators.h"


namespace petibm
{
namespace applications
{

/**
 * \class NavierStokesSolver
 * \brief Solve the incompressible Navier-Stokes equations with a projection 
 *        method Perot (1993).
 */
class NavierStokesSolver
{
public:
	/**
	 * \brief Default constructor.
	 */
	NavierStokesSolver();
	
	/**
	 * \brief Constructor. Set references to the mesh, flow conditions,
	 *        and simulation parameters.
	 *
	 * \param mesh Structured Cartesian mesh
	 * \param flow Flow conditions
	 * \param parameters Simulation parameters
	 */
	NavierStokesSolver(utilities::CartesianMesh &mesh,
	                   utilities::FlowDescription &flow,
	                   utilities::SimulationParameters &parameters);

	/**
	 * \brief Default destructor.
	 */
	~NavierStokesSolver();

	/**
	 * \brief Initialize vectors, operators, and linear solvers.
	 */
	PetscErrorCode initialize();

	/**
	 * \brief Advance in time.
	 */
	PetscErrorCode solve();

	/**
	 * \brief Write the solution into a file.
	 *
	 * \param directory Directory where to save the file
	 * \param fileName Name of the file to save (without the extension)
	 */
	PetscErrorCode write(std::string directory, std::string fileName);

	/**
	 * \brief Write number of iterations executed by each solver at current time
	 *        step.
	 *
	 * \param timeIndex Time-step index
	 * \param filePath Path of the file to write in
	 */
	PetscErrorCode writeIterations(int timeIndex, std::string filePath);

	/**
	 * \brief Destroy PETSc objects (vectors and matrices) and linear solvers.
	 */
	PetscErrorCode finalize();

private:
	utilities::CartesianMesh mesh;  ///< Structured Cartesian mesh
	utilities::FlowDescription flow;  ///< Flow conditions
	utilities::SimulationParameters parameters;  ///< Simulation parameters
	utilities::Solutions solution;  ///< Velocity and pressure fields
	utilities::Boundary bc;  ///< Information about the domain boundaries
	utilities::TimeIntegration convection,  ///< Time scheme for the convective terms
	                           diffusion;  ///< Time scheme for the diffusive terms

	Mat L,  ///< Laplacian operator
	    LCorrection,  ///< Laplacian correction for boundary conditions
	    G,  ///< Gradient operator
	    D,  ///< Divergence operator
	    DCorrection,  ///< Divergence correction for boundary conditions
	    N,  ///< Linear convection operator
	    I,  ///< Identity matrix
	    R,
	    RInv,
	    M,
	    MHat,
	    BNHat,
	    A,  ///< Matrix resulting from implicit treatment
	    BNG, ///< Projection operator
	    DBNG;  ///< Poisson matrix

	Vec phi,  ///< Pressure-correction vector
	    rhs1,  ///< Right-hand side vector of the velocity system
	    bc1,  ///< Boundary terms for the velocity system
	    rhs2,  ///< Right-hand side vector of the Poisson system
	    gradP;  ///< Pressure-gradient vector
	std::vector<Vec> Conv,  ///< Convective terms from previous time steps
	                 Diff;  ///< Diffusive terms from previous time steps

	std::shared_ptr<linsolvers::LinSolver> vSolver,  ///< Velocity linear solver
	                                       pSolver;  ///< Poisson linear solver

	PetscLogStage stageInitialize,  ///< Log initialize phase
	              stageRHSVelocity,  ///< Log RHS of velocity system
	              stageSolveVelocity,  ///< Log velocity solve
	              stageRHSPoisson,  ///< Log RHS of Poisson system
	              stageSolvePoisson,  ///< Log Poisson solve
	              stageProjectionStep,  ///< Log projection step
	              stageWrite;  ///< Log write phase

	/**
	 * \brief Assembles the different operators and matrices.
	 */
	PetscErrorCode assembleOperators();

	/**
	 * \brief Assemble the RHS vector of the velocity system.
	 */
	PetscErrorCode assembleRHSVelocity();

	/**
	 * \brief Solve the velocity system.
	 */
	PetscErrorCode solveVelocity();

	/**
	 * \brief Assemble the RHS vector of the Poisson system.
	 */
	PetscErrorCode assembleRHSPoisson();

	/**
	 * \brief Solve the Poisson system.
	 */
	PetscErrorCode solvePoisson();

	/**
	 * \brief Project the velocity field onto the divergence-free space and update
	 *        the pressure field.
	 */
	PetscErrorCode projectionStep();
};

} // end of namespace applications
} // end of namespace petibm
