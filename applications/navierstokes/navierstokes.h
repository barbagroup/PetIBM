/**
 * \file navierstokes.h
 * \brief Definition of the class \c NavierStokesSolver.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 * \see nssolver
 * \ingroup nssolver
 */

#pragma once

#include <yaml-cpp/yaml.h>

#include <petibm/boundary.h>
#include <petibm/linsolver.h>
#include <petibm/mesh.h>
#include <petibm/operators.h>
#include <petibm/probes.h>
#include <petibm/solution.h>
#include <petibm/timeintegration.h>

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
    NavierStokesSolver(const petibm::type::Mesh &mesh,
                       const petibm::type::Boundary &bc,
                       const YAML::Node &node);

    /**
     * \brief Default destructor.
     */
    ~NavierStokesSolver();

    /**
     * \brief Manually destroy data.
     */
    PetscErrorCode destroy();

    /**
     * \brief Initialize vectors, operators, and linear solvers.
     *
     * \param mesh [in] Structured Cartesian mesh.
     * \param bc [in] Data structure with the boundary conditions.
     * \param node [in] YAML node with configuration parameters.
     */
    PetscErrorCode initialize(const petibm::type::Mesh &mesh,
                              const petibm::type::Boundary &bc,
                              const YAML::Node &node);

    /**
     * \brief Advance the solution in time by one time step.
     */
    PetscErrorCode advance();

    /**
     * \brief Write the solution fields into a file.
     *
     * \param t [in] Time.
     * \param filePath [in] Path of the file to write into.
     */
    PetscErrorCode write(const PetscReal &t, const std::string &filePath);

    /**
     * \brief Write into a file all data required to restart a simulation.
     *
     * If the directory already contains the solution fields at the present
     * time, the method only writes the convective and diffusive terms into it.
     *
     * \param t [in] Time.
     * \param filePath [in] Path of the file to write into.
     */
    PetscErrorCode writeRestartData(const PetscReal &t,
                                    const std::string &filePath);

    /**
     * \brief Read from a file all data required to restart a simulation.
     *
     * \param filePath [in] Path of the file to read from.
     * \param t [out] Time.
     */
    PetscErrorCode readRestartData(const std::string &filePath, PetscReal &t);

    /**
     * \brief Initialize PetscViewer objects for ASCII files.
     *
     * The method can be used to initialize a viewer for the iterations file or
     * the PETSc Log summary file.
     *
     * \param filePath [in] Path of the file to write into.
     * \param mode [in] Either FILE_MODE_WRITE (default) or FILE_MODE_APPEND.
     */
    PetscErrorCode initializeASCIIFiles(
        const std::string &filePath,
        const PetscFileMode &mode = FILE_MODE_WRITE);

    /**
     * \brief Output the number of iterations for each solver at current time
     * step.
     *
     * Requires initializing a PetscViewer object with the method
     * `initializeASCIIFiles` for the given `filePath`.
     *
     * \param timeIndex [in] Time-step index.
     * \param filePath [in] Path of the file to write in.
     */
    PetscErrorCode writeIterations(const int &timeIndex,
                                   const std::string &filePath);

    /**
     * \brief Write the time value into HDF5 file.
     *
     * The time value is written as an attribute of the pressure (dataset `p`),
     * which must be existing.
     *
     * \param t [in] Time.
     * \param filePath [in] Path of the file to write in.
     */
    PetscErrorCode writeTimeHDF5(const PetscReal &t,
                                 const std::string &filePath);

    /**
     * \brief Read the time value from HDF5 file.
     *
     * The time value is an attribute of the pressure (dataset `p`),
     * which must be existing.
     *
     * \param filePath [in] Path of the file to read from.
     * \param t [out] Time.
     */
    PetscErrorCode readTimeHDF5(const std::string &filePath, PetscReal &t);

    /**
     * \brief Monitor the solutions at the probes.
     *
     * \param t [in] Time
     * \param ite [in] Time-step index
     * \return PetscErrorCode
     */
    PetscErrorCode monitorProbes(const PetscReal &t, const PetscInt &ite);

protected:
    /** \brief A reference to the YAML::Node passed in. */
    YAML::Node settings;

    /** \brief Structured Cartesian mesh. */
    petibm::type::Mesh mesh;

    /** \brief Information about the domain boundaries. */
    petibm::type::Boundary bc;

    /** \brief Velocity and pressure fields. */
    petibm::type::Solution solution;

    /** \brief Time scheme for the convective terms. */
    petibm::type::TimeIntegration convCoeffs;

    /** \brief Time scheme for the diffusive terms. */
    petibm::type::TimeIntegration diffCoeffs;

    /** \brief Velocity linear solver. */
    petibm::type::LinSolver vSolver;

    /** \brief Poisson linear solver. */
    petibm::type::LinSolver pSolver;

    /** \brief Probes to monitor the solution. */
    std::vector<petibm::type::Probe> probes;

    /** \brief A copy of time-step size. */
    PetscReal dt;

    /** \brief A copy of viscosity. */
    PetscReal nu;

    /** \brief Laplacian operator. */
    Mat L;

    /** \brief Laplacian correction for boundary conditions. */
    Mat LCorrection;

    /** \brief Gradient operator. */
    Mat G;

    /** \brief Divergence operator. */
    Mat D;

    /** \brief Divergence correction for boundary conditions. */
    Mat DCorrection;

    /** \brief Non-linear convection operator (matrix-free). */
    Mat N;

    /** \brief Matrix resulting from implicit treatment. */
    Mat A;

    /** \brief Projection operator. */
    Mat BNG;

    /** \brief Poisson matrix. */
    Mat DBNG;

    /** \brief Pressure-correction vector. */
    Vec dP;

    /** \brief Boundary terms for the velocity system. */
    Vec bc1;

    /** \brief Right-hand side vector of the velocity system. */
    Vec rhs1;

    /** \brief Right-hand side vector of the Poisson system. */
    Vec rhs2;

    /** \brief Convective terms from previous time steps. */
    std::vector<Vec> conv;

    /** \brief Diffusive terms from previous time steps. */
    std::vector<Vec> diff;

    /** \brief A Bool indicating if we'll pin a reference pressure point. */
    PetscBool isRefP;

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

    /** \brief Log monitor phase. */
    PetscLogStage stageMonitor;

    /** \brief A dictionary mapping file path to PetscViewer objects. */
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
     *         update the pressure field.
     */
    virtual PetscErrorCode projectionStep();

    /** \brief Create operators.  */
    virtual PetscErrorCode createOperators();

    /** \brief Create vectors.  */
    virtual PetscErrorCode createVectors();

    /** \brief Set null space or apply reference point.  */
    virtual PetscErrorCode setNullSpace();

};  // NavierStokesSolver
