/**
 * \file navierstokes.h
 * \brief Definition of the class \c NavierStokesSolver.
 * \see nssolver
 * \ingroup nssolver
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
 *        method (Perot 1993).
 * \see nssolver
 * \ingroup nssolver
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
    ~NavierStokesSolver();


    /**
     * \brief Manually destroy data.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode destroy();

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
     * be written in. Otherwise, solutions and extra data will all be written in.
     *
     * \param filePath [in] path of the file to save (without the extension)
     */
    PetscErrorCode writeRestartData(const std::string &filePath);
    
    /**
     * \brief Read data that are required for restarting sessions.
     * 
     * \param filePath [in] path of the file to save (without the extension)
     */
    PetscErrorCode readRestartData(const std::string &filePath);

    /**
     * \brief Initialize viewers for ASCII files, such as iteration log.
     *
     * \param filePath [in] a string indicating the path to the file.
     * \param mode [in] either FILE_MODE_WRITE (default) or FILE_MODE_APPEND.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode initializeASCIIFiles(const std::string &filePath,
            const PetscFileMode &mode=FILE_MODE_WRITE);

    /**
     * \brief Write number of iterations executed by each solver at current time
     *        step.
     *
     * For a given `filePath`, an `initializeASCIIFiles` should be called for 
     * this `filePath` prior any call to `writeIterations`.
     *
     * \param timeIndex [in] Time-step index
     * \param filePath [in] Path of the file to write in
     */
    PetscErrorCode writeIterations(const int &timeIndex, const std::string &filePath);

protected:
    
    
    /** \brief A reference to the YAML::Node passed in. */
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
    

    
    /** \brief A copy of time-step size. */
    PetscReal   dt;
    
    /** \brief A copy of viscosity. */
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
    

    
    /** \brief A Bool indicating if we'll pin a reference pressure point. */
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



    /** \brief A dictionary mapping file path to PetscViewers. */
    std::map<std::string, PetscViewer> asciiViewers;
    

    

    /** \brief Assemble the RHS vector of the velocity system.  */
    virtual PetscErrorCode assembleRHSVelocity();

    /** \brief Solve the velocity system.  */
    virtual PetscErrorCode solveVelocity();

    /** \brief Assemble the RHS vector of the Poisson system.  */
    virtual PetscErrorCode assembleRHSPoisson();

    /** \brief Solve the Poisson system.  */
    virtual PetscErrorCode solvePoisson();

    /** \brief Project the velocity field onto the divergence-free space and 
     *         update the pressure field.  */
    virtual PetscErrorCode projectionStep();
    
    /** \brief Create operators.  */
    virtual PetscErrorCode createOperators();
    
    /** \brief Create vectors.  */
    virtual PetscErrorCode createVectors();
    
    /** \brief Set null space or apply reference point.  */
    virtual PetscErrorCode setNullSpace();
}; // NavierStokesSolver
