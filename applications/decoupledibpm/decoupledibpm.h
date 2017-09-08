/**
 * \file decoupledibpm.h
 * \brief Definition of the class \c DecoupledIBPMSolver.
 */

#pragma once

#include "petibm/cartesianmesh.h"
#include "petibm/flowdescription.h"
#include "petibm/simulationparameters.h"
#include "petibm/bodypack.h"
#include "petibm/solutions.h"
#include "petibm/boundary.h"
#include "petibm/timeintegration.h"
#include "petibm/linsolver.h"
#include "petibm/operators.h"

using namespace petibm;


/**
 * \class DecoupledIBPMSolver
 * \brief Solve the incompressible Navier-Stokes equations with a decoupled
 *        version of the Immersed-Boundary Projection Method (Li et al., 2016).
 */
class DecoupledIBPMSolver
{
public:
	/**
	 * \brief Default constructor.
	 */
	DecoupledIBPMSolver();
	
	/**
	 * \brief Constructor. Set references to the mesh, flow conditions,
	 *        and simulation parameters.
	 *
	 * \param mesh Structured Cartesian mesh
	 * \param flow Flow conditions
	 * \param parameters Simulation parameters
	 */
	DecoupledIBPMSolver(utilities::CartesianMesh &mesh,
	                    utilities::FlowDescription &flow,
	                    utilities::SimulationParameters &parameters,
	                    utilities::BodyPack &bodies);

	/**
	 * \brief Default destructor.
	 */
	~DecoupledIBPMSolver();

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
	 * \brief Write the integrated forces acting on the bodies into a file.
	 *
	 * \param time Time value
	 * \param directory Directory where to save the file
	 * \param fileName Name of the file to save (without the extension)
	 */
	PetscErrorCode writeIntegratedForces(
			int time, std::string directory, std::string fileName);

	/**
	 * \brief Destroy PETSc objects (vectors and matrices) and linear solvers.
	 */
	PetscErrorCode finalize();

private:
	utilities::CartesianMesh mesh;  ///< Structured Cartesian mesh
	utilities::FlowDescription flow;  ///< Flow conditions
	utilities::SimulationParameters parameters;  ///< Simulation parameters
	utilities::BodyPack bodies;
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
	    EHat,
	    HHat,
	    BNHHat,
	    EBNHHat,
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
	Vec f,
	    df,
	    Eu,
	    Hf;

	std::shared_ptr<linsolvers::LinSolver> vSolver,  ///< Velocity linear solver
	                                       pSolver,  ///< Poisson linear solver
	                                       fSolver;  ///< Forces linear solver

	PetscLogStage stageInitialize,  ///< Log initialize phase
	              stageRHSVelocity,  ///< Log RHS of velocity system
	              stageSolveVelocity,  ///< Log velocity solve
	              stageRHSPoisson,  ///< Log RHS of Poisson system
	              stageSolvePoisson,  ///< Log Poisson solve
	              stageRHSForces,  ///< Log RHS of forces system
	              stageSolveForces,  ///< Log forces solver
	              stageProjectionStep,  ///< Log projection step
	              stageIntegrateForces,  ///< Log force integration
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
	 * \brief Assemble the RHS vector of the system for the boundary forces.
	 */
	PetscErrorCode assembleRHSForces();

	/**
	 * \brief Solve the system for the boundary forces.
	 */
	PetscErrorCode solveForces();

	/**
	 * \brief Project the velocity field onto the divergence-free space and update
	 *        the pressure field.
	 */
	PetscErrorCode projectionStep();

}; // DecoupledIBPMSolver
