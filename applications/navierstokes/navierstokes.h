/**
 * \file navierstokes.h
 * \brief Definition of the class \c NavierStokesSolver.
 */

#pragma once

// YAML-CPP
# include <yaml-cpp/yaml.h>

// PetIBM
# include <petibm/mesh.h>
# include <petibm/solution.h>
# include <petibm/boundary.h>
# include <petibm/timeintegration.h>
# include <petibm/linsolver.h>
# include <petibm/operators.h>


/**
 * \class NavierStokesSolver
 * \brief Solve the incompressible Navier-Stokes equations with a projection 
 *        method Perot (1993).
 */
class NavierStokesSolver
{
public:
    
    /** \brief Default constructor.  */
    NavierStokesSolver() = default;

    /**
     * \brief Constructor; Set references to the mesh and boundary conditions.
     *
     * \param mesh [in] a type::Mesh object.
     * \param bc [in] a type::Boundary object.
     * \param node [in] YAML::Node containing settings.
     */
    NavierStokesSolver(
            const petibm::type::Mesh &mesh,
            const petibm::type::Boundary &bc,
            const YAML::Node &node);

    /**
     * \brief Default destructor.
     */
    ~NavierStokesSolver() = default;

    /**
     * \brief Initialize vectors, operators, and linear solvers.
     */
    PetscErrorCode initialize(
            const petibm::type::Mesh &mesh,
            const petibm::type::Boundary &bc,
            const YAML::Node &node);

    /**
     * \brief Advance in time for one step.
     */
    PetscErrorCode advance();

    /**
     * \brief Write the solution into a file.
     *
     * \param filePath [in] path of the file to save (without the extension)
     */
    PetscErrorCode write(const std::string &filePath);
    
    /**
     * \brief Write the extra data that are required for restarting sessions.
     * 
     * If the file already has solutions in it, only extra necessary data will
     * be writen in. Otherwise, solutions and extra data will all be writen in.
     *
     * \param filePath [in] path of the file to save (without the extension)
     */
    PetscErrorCode writeRestartData(const std::string &filePath);
    
    /**
     * \brief read data that are required for restarting sessions.
     * 
     * \param filePath [in] path of the file to save (without the extension)
     */
    PetscErrorCode readRestartData(const std::string &filePath);

    /**
     * \brief Write number of iterations executed by each solver at current time
     *        step.
     *
     * \param timeIndex Time-step index
     * \param filePath Path of the file to write in
     */
    PetscErrorCode writeIterations(
            const int &timeIndex, const std::string &filePath);

    /**
     * \brief Destroy PETSc objects (vectors and matrices) and linear solvers.
     */
    PetscErrorCode finalize();

    
protected:
    
    
    /** \brief a reference to the YAML::Node passed in. */
    YAML::Node                      settings;
    
    /** \brief Structured Cartesian mesh. */
    petibm::type::Mesh              mesh; 
    
    /** \brief Information about the domain boundaries. */
    petibm::type::Boundary          bc; 
    
    /** \brief Velocity and pressure fields. */
    petibm::type::Solution          solution; 
    
    /** \brief Time scheme for the convective terms. */
    petibm::type::TimeIntegration   convCoeffs; 
    
    /** \brief Time scheme for the diffusive terms. */
    petibm::type::TimeIntegration   diffCoeffs; 

    /** \brief Velocity linear solver. */
    petibm::type::LinSolver         vSolver; 
    
    /** \brief Poisson linear solver. */
    petibm::type::LinSolver         pSolver; 
    

    
    /** \brief a copy of time-step size. */
    PetscReal   dt;
    
    /** \brief a copy of viscosity. */
    PetscReal   nu;
    

    
    /** \brief Laplacian operator. */
    Mat         L; 

    /** \brief Laplacian correction for boundary conditions. */
    Mat         LCorrection; 

    /** \brief Gradient operator. */
    Mat         G; 

    /** \brief Divergence operator. */
    Mat         D; 

    /** \brief Divergence correction for boundary conditions. */
    Mat         DCorrection; 

    /** \brief Non-linear convection operator (matrix-free). */
    Mat         N; 

    /** \brief Matrix resulting from implicit treatment. */
    Mat         A; 

    /** \brief Projection operator. */
    Mat         BNG;

    /** \brief Poisson matrix. */
    Mat         DBNG; 
    
    

    /** \brief Pressure-correction vector. */
    Vec         dP; 
    
    /** \brief Boundary terms for the velocity system. */
    Vec         bc1; 
    
    /** \brief Right-hand side vector of the velocity system. */
    Vec         rhs1; 
    
    /** \brief Right-hand side vector of the Poisson system. */
    Vec         rhs2; 
    
    /** \brief Convective terms from previous time steps. */
    std::vector<Vec> conv; 
    
    /** \brief Diffusive terms from previous time steps. */
    std::vector<Vec> diff; 
    

    
    /** \brief a bool indicating if we'll pin a reference pressure point. */
    PetscBool   isRefP;
    
    

    /** \brief Log initialize phase. */
    PetscLogStage stageInitialize;
    
    /** \brief Log RHS of velocity system. */
    PetscLogStage stageRHSVelocity;
    
    /** \brief Log velocity solve. */
    PetscLogStage stageSolveVelocity;
    
    /** \brief Log RHS of Poisson system. */
    PetscLogStage stageRHSPoisson;
    
    /** \brief Log Poisson solve. */
    PetscLogStage stageSolvePoisson;
    
    /** \brief Log projection step. */
    PetscLogStage stageProjectionStep;
    
    /** \brief Log write phase. */
    PetscLogStage stageWrite;
    

    

    /** \brief Assemble the RHS vector of the velocity system.  */
    PetscErrorCode assembleRHSVelocity();

    /** \brief Solve the velocity system.  */
    PetscErrorCode solveVelocity();

    /** \brief Assemble the RHS vector of the Poisson system.  */
    PetscErrorCode assembleRHSPoisson();

    /** \brief Solve the Poisson system.  */
    PetscErrorCode solvePoisson();

    /** \brief Project the velocity field onto the divergence-free space and 
     *         update the pressure field.  */
    PetscErrorCode projectionStep();
    
    
    PetscErrorCode createOperators();
    PetscErrorCode createVectors();
    PetscErrorCode setNullSpace();
};
