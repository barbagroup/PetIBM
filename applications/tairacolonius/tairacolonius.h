/**
 * \file tairacolonius.h
 * \brief Definition of the class \c TairaColoniusSolver.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 * \see tairacolonius
 * \ingroup tairacolonius
 */

#pragma once

// Navier-Stokes solver
# include "../navierstokes/navierstokes.h"

// PetIBM
# include <petibm/bodypack.h>


/**
 * \class TairaColoniusSolver
 * \brief Immersed-boundary method proposed by Taira and Colonius (2007).
 * \see tairacolonius, NavierStokesSolver
 * \ingroup tairacolonius
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

    /** \brief Default destructor.  */
    ~TairaColoniusSolver();

    /** \brief Manually destroy data.  */
    PetscErrorCode destroy();

    /** \brief Initialize vectors, operators, and linear solvers.  */
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
     * be written in. Otherwise, solutions and extra data will all be written in.
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
     * \param t [in] Time value
     * \param filePath [in] Name of the file to save.
     */
    PetscErrorCode writeIntegratedForces(
            const PetscReal &t, const std::string &filePath);

    
protected:
    
    /** \brief A reference to immersed bodies. */
    petibm::type::BodyPack      bodies;
    
    /** \brief Combination of pressure and forces. */
    Vec P;
    
    /** \brief PETSc IS objects indicating which entries in phi belonging to 
     *         pressure or forces. */
    IS isDE[2];
    
    /** \brief Log force integration. */
    PetscLogStage stageIntegrateForces; 
    


    /** \brief Assemble the RHS vector of the Poisson system.  */
    virtual PetscErrorCode assembleRHSPoisson();
    
    /** \brief Create operators.  */
    virtual PetscErrorCode createOperators();
    
    /** \brief Create vectors.  */
    virtual PetscErrorCode createVectors();
    
    /** \brief Set null space or apply reference point.  */
    virtual PetscErrorCode setNullSpace();
}; // TairaColoniusSolver
