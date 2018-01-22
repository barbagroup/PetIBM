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

    // public methods that don't change
    using NavierStokesSolver::advance;
    using NavierStokesSolver::initializeASCIIFiles;
    using NavierStokesSolver::writeIterations;
    using NavierStokesSolver::writeTimeHDF5;
    using NavierStokesSolver::readTimeHDF5;
    
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
    ~TairaColoniusSolver();

    /** \brief manually destroy data.  */
    PetscErrorCode destroy();

    /**
     * \brief Initialize vectors, operators, and linear solvers.
     */
    PetscErrorCode initialize(
            const petibm::type::Mesh &mesh,
            const petibm::type::Boundary &bc,
            const petibm::type::BodyPack &bodies,
            const YAML::Node &node);

    /**
     * \brief Write the solution into a file.
     *
     * \param t [in] time
     * \param filePath [in] path of the file to save (without the extension)
     */
    PetscErrorCode write(const PetscReal &t, const std::string &filePath);
    
    /**
     * \brief Write the extra data that are required for restarting sessions.
     * 
     * If the file already has solutions in it, only extra necessary data will
     * be writen in. Otherwise, solutions and extra data will all be writen in.
     *
     * \param t [in] time
     * \param filePath [in] path of the file to save (without the extension)
     */
    PetscErrorCode writeRestartData(
      const PetscReal &t, const std::string &filePath);
    
    /**
     * \brief read data that are required for restarting sessions.
     * 
     * \param filePath [in] path of the file to save (without the extension)
     * \param t [out] time
     */
    PetscErrorCode readRestartData(const std::string &filePath, PetscReal &t);

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
    
    /** \brief combination of pressure and forces. */
    Vec P;
    
    /** \brief PETSc IS objects indicating which entries in phi belonging to 
     *         pressure or forces. */
    IS isDE[2];
    
    /** \brief Log force integration. */
    PetscLogStage stageIntegrateForces; 
    


    /** \brief Assemble the RHS vector of the Poisson system.  */
    virtual PetscErrorCode assembleRHSPoisson();
    
    /** \brief create operators.  */
    virtual PetscErrorCode createOperators();
    
    /** \brief create vectors.  */
    virtual PetscErrorCode createVectors();
    
    /** \brief set null space or apply reference point.  */
    virtual PetscErrorCode setNullSpace();
};
