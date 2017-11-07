/**
 * \file navierstokes.h
 * \brief Definition of the class \c TairaColoniusSolver.
 */

#pragma once

// Navier-Stokes solver
# include "../navierstokes/navierstokes.h"

// PetIBM
# include <petibm/bodypack.h>


/**
 * \class TairaColoniusSolver
 * \brief Taira and Colonius (2007).
 */
class TairaColoniusSolver : protected NavierStokesSolver
{
public:
    
    /** \brief Default constructor.  */
    TairaColoniusSolver() = default;

    /**
     * \brief Constructor; Set references to the mesh, boundary conditions, and
     *        immersed bodies.
     *
     * \param mesh [in] a type::Mesh object.
     * \param bc [in] a type::Boundary object.
     * \param bodies [in] a type::BodyPack object.
     * \param node [in] YAML::Node containing settings.
     */
    TairaColoniusSolver(
            const petibm::type::Mesh &mesh,
            const petibm::type::Boundary &bc,
            const petibm::type::BodyPack &bodies,
            const YAML::Node &node);

    /**
     * \brief Default destructor.
     */
    ~TairaColoniusSolver() = default;

    /**
     * \brief Initialize vectors, operators, and linear solvers.
     */
    PetscErrorCode initialize(
            const petibm::type::Mesh &mesh,
            const petibm::type::Boundary &bc,
            const petibm::type::BodyPack &bodies,
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
     * \brief Write the integrated forces acting on the bodies into a ASCII file.
     *
     * \param time [in] Time value
     * \param fileName [in] Name of the file to save.
     */
    PetscErrorCode writeIntegratedForces(
            const PetscReal &t, const std::string &filePath);

    /**
     * \brief Destroy PETSc objects (vectors and matrices) and linear solvers.
     */
    PetscErrorCode finalize();

    
protected:
    
    /** \brief a reference to immersed bodies. */
    petibm::type::BodyPack      bodies;
    
    
    
    /** \brief a bool indicating if we'll pin a reference pressure point. */
    PetscBool   isRefP;
    
    
    /** \brief solution of Lagragian force at time-step n. */
    Vec f;
    
    /** \brief combination of pressure and forces. */
    Vec phi;
    
    /** \brief increment of phi. */
    Vec dphi;
    
    /** \brief boundary correction from divergence operator. */
    Vec bc2;
    

    IS isDE[2];
    
    
    
    /** \brief Log force integration. */
    PetscLogStage stageIntegrateForces; 
    
    

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
