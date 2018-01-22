/**
 * \file decoupledibpm.h
 * \brief Definition of the class \c DecoupledIBPMSolver.
 */

#pragma once

// Navier-Stokes solver
# include "../navierstokes/navierstokes.h"

// additional PetIBM headers
# include <petibm/bodypack.h>


/**
 * \class DecoupledIBPMSolver
 * \brief Solve the incompressible Navier-Stokes equations with a decoupled
 *        version of the Immersed-Boundary Projection Method (Li et al., 2016).
 */
class DecoupledIBPMSolver : protected NavierStokesSolver
{
public:

    // public members that don't change
    using NavierStokesSolver::write;
    using NavierStokesSolver::initializeASCIIFiles;
    using NavierStokesSolver::writeTimeHDF5;
    using NavierStokesSolver::readTimeHDF5;
    
    /** \brief Default constructor. */
    DecoupledIBPMSolver() = default;
    
    /**
     * \brief Constructor; Set references to the mesh, boundary conditions, and
     *        immersed bodies.
     *
     * \param mesh [in] a type::Mesh object.
     * \param bc [in] a type::Boundary object.
     * \param bodies [in] a type::BodyPack object.
     * \param node [in] YAML::Node containing settings.
     */
    DecoupledIBPMSolver(
            const petibm::type::Mesh &mesh,
            const petibm::type::Boundary &bc,
            const petibm::type::BodyPack &bodies,
            const YAML::Node &node);

    /** \brief Default destructor. */
    ~DecoupledIBPMSolver();

    /** \brief manually destroy data. */
    PetscErrorCode destroy();

    /** \brief Initialize vectors, operators, and linear solvers. */
    PetscErrorCode initialize(
            const petibm::type::Mesh &mesh,
            const petibm::type::Boundary &bc,
            const petibm::type::BodyPack &bodies,
            const YAML::Node &node);

    /** \brief Advance in time. */
    PetscErrorCode advance();
    
    /**
     * \brief Write the extra data that are required for restarting sessions.
     * 
     * If file already exists, only extra necessary data will
     * be writen in. Otherwise, solutions and extra data will all be writen in.
     *
     * \param t [in] time
     * \param filePath [in] path of the file to save (without extension)
     */
    PetscErrorCode writeRestartData(
      const PetscReal &t, const std::string &filePath);
    
    /**
     * \brief read data that are required for restarting sessions.
     * 
     * \param filePath [in] path of the file to save (without extension)
     * \param t [out] time
     */
    PetscErrorCode readRestartData(const std::string &filePath, PetscReal &t);

    /**
     * \brief Write number of iterations executed by each solver at current time
     *        step (to an ASCII file).
     *
     * \param timeIndex [in] Time-step index
     * \param filePath [in] Path of the file to write in
     */
    PetscErrorCode writeIterations(
            const int &timeIndex, const std::string &filePath);

    /**
     * \brief Write the integrated forces acting on the bodies into a ASCII file.
     *
     * \param time [in] Time value
     * \param fileName [in] Name of the file to save.
     */
    PetscErrorCode writeIntegratedForces(
            const PetscReal &t, const std::string &filePath);

protected:
    
    /** \brief a reference to immersed bodies. */
    petibm::type::BodyPack      bodies;
    
    /** \brief linear solver object for force solver. */
    petibm::type::LinSolver     fSolver;
    

    /** \brief operator interpolating Lagrangian forces to Eulerian forces. */
    Mat H;
    
    /** \brief operator interpolating Eulerian forces to Lagrangian forces. */
    Mat E;
    
    /** \brief coefficient matrix of the force system. */
    Mat EBNH;
    
    /** \brief operator projecting force to intermediate velocity field. */
    Mat BNH;
    

    /** \brief right-hand-side of force system. */
    Vec Eu;
    
    /** \brief solution of Lagragian force at time-step n. */
    Vec f;
    
    /** \brief increment of force from time-step n to n+1. */
    Vec df;


    /** \brief Log RHS of forces system. */
    PetscLogStage stageRHSForces; 
    
    /** \brief Log forces solver. */
    PetscLogStage stageSolveForces; 
    
    /** \brief Log force integration. */
    PetscLogStage stageIntegrateForces; 
    

    /** \brief Assemble the RHS vector of the velocity system.  */
    virtual PetscErrorCode assembleRHSVelocity();

    /** \brief Assemble the RHS vector of the Poisson system. */
    virtual PetscErrorCode assembleRHSPoisson();

    /** \brief Assemble the RHS vector of the system for the boundary forces. */
    virtual PetscErrorCode assembleRHSForces();

    /** \brief Solve the system for the boundary forces. */
    virtual PetscErrorCode solveForces();

    /** \brief Project the velocity to divergence-free space, update
     *         pressure field, and update force.  */
    virtual PetscErrorCode projectionStep();

    /** \brief Assembles operators and matrices. */
    PetscErrorCode createExtraOperators();

    /** \brief create vectors. */
    PetscErrorCode createExtraVectors();

}; // DecoupledIBPMSolver
