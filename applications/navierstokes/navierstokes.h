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
    /** \brief Default constructor. */
    NavierStokesSolver() = default;

    /** \brief Constructor. Initialize the Navier-Stokes solver.
     *
     * \param world [in] MPI communicator
     * \param node [in] YAML configuration settings
     */
    NavierStokesSolver(const MPI_Comm &world,
                       const YAML::Node &node);

    /** \brief Default destructor. */
    ~NavierStokesSolver();

    /** \brief Manually destroy data. */
    PetscErrorCode destroy();

    /** \brief Initialize the Navier-Stokes solver.
     *
     * Create a structured Cartesian mesh and write it into a HDF5 file.
     * Create the initial field solution.
     * Create the linear solvers (along with their vectors and operators).
     * Create probes to monitor the solution in some user-defined regions.
     * Record some simulation parameters.
     *
     * \param world [in] MPI communicator
     * \param node [in] YAML configuration settings
     * \return PetscErrorCode
     */
    PetscErrorCode init(const MPI_Comm &world, const YAML::Node &node);

    /** \brief Read or write initial data. */
    PetscErrorCode ioInitialData();

    /** \brief Advance the solution by one time step. */
    PetscErrorCode advance();

    /** \brief Write solution and solver info to files. */
    PetscErrorCode write();

    /** \brief Evaluate if the simulation is finished. */
    bool finished();

protected:
    /** \brief MPI communicator. */
    MPI_Comm comm;

    /** \brief Size of the MPI communicator. */
    PetscMPIInt commSize;

    /** \brief Rank of the process in the MPI communicator. */
    PetscMPIInt commRank;

    /** \brief YAML configuration settings. */
    YAML::Node config;

    /** \brief Structured Cartesian mesh object. */
    petibm::type::Mesh mesh;

    /** \brief Information about the domain boundaries. */
    petibm::type::Boundary bc;

    /** \brief Data object holding the velocity and pressure fields. */
    petibm::type::Solution solution;

    /** \brief Time scheme for the convective terms. */
    petibm::type::TimeIntegration convCoeffs;

    /** \brief Time scheme for the diffusion terms. */
    petibm::type::TimeIntegration diffCoeffs;

    /** \brief Velocity linear solver. */
    petibm::type::LinSolver vSolver;

    /** \brief Poisson linear solver. */
    petibm::type::LinSolver pSolver;

    /** \brief Probes to monitor the solution. */
    std::vector<petibm::type::Probe> probes;

    /** \brief Time-step size. */
    PetscReal dt;

    /** \brief Time-step index. */
    PetscInt ite;

    /** \brief Time value. */
    PetscReal t;

    /** \brief Starting time-step index. */
    PetscInt nstart;

    /** \brief Number of time steps to compute. */
    PetscInt nt;

    /** \brief Frequency at which the solution fields are written to files. */
    PetscInt nsave;

    /** \brief Frequency at which data to restart are written to files. */
    PetscInt nrestart;

    /** \brief Viscous diffusion coefficient. */
    PetscReal nu;

    /** \brief Laplacian operator. */
    Mat L;

    /** \brief Laplacian correction operator for boundary conditions. */
    Mat LCorrection;

    /** \brief Gradient operator. */
    Mat G;

    /** \brief Divergence operator. */
    Mat D;

    /** \brief Divergence correction for boundary conditions. */
    Mat DCorrection;

    /** \brief Convective operator (matrix-free). */
    Mat N;

    /** \brief Implicit operator for the velocity solver. */
    Mat A;

    /** \brief Projection operator. */
    Mat BNG;

    /** \brief Poisson operator. */
    Mat DBNG;

    /** \brief Pressure-correction vector. */
    Vec dP;

    /** \brief Inhomogeneous boundary terms for the velocity system. */
    Vec bc1;

    /** \brief Right-hand side vector of the velocity system. */
    Vec rhs1;

    /** \brief Right-hand side vector of the Poisson system. */
    Vec rhs2;

    /** \brief Explicit convective terms. */
    std::vector<Vec> conv;

    /** \brief Explicit diffusion terms. */
    std::vector<Vec> diff;

    /** \brief True if we pin the pressure at a reference point. */
    PetscBool isRefP;

    /** \brief Log stage for the initialization phase. */
    PetscLogStage stageInitialize;

    /** \brief Log stage for assembling the RHS of the velocity system. */
    PetscLogStage stageRHSVelocity;

    /** \brief Log stage for solving the velocity system. */
    PetscLogStage stageSolveVelocity;

    /** \brief Log stage for assembling the RHS of the Poisson system. */
    PetscLogStage stageRHSPoisson;

    /** \brief Log stage for solving the Poisson system. */
    PetscLogStage stageSolvePoisson;

    /** \brief Log stage for projecting the velocity field. */
    PetscLogStage stageProjectionStep;

    /** \brief Log stage when write the solution fields. */
    PetscLogStage stageWrite;

    /** \brief Log stage for monitoring user-defined regions of the domain. */
    PetscLogStage stageMonitor;

    /** \brief ASCII PetscViewer object to output solvers info. */
    PetscViewer solversViewer;

    /** \brief Assemble the RHS vector of the velocity system. */
    virtual PetscErrorCode assembleRHSVelocity();

    /** \brief Solve the velocity system. */
    virtual PetscErrorCode solveVelocity();

    /** \brief Assemble the RHS vector of the Poisson system. */
    virtual PetscErrorCode assembleRHSPoisson();

    /** \brief Solve the Poisson system. */
    virtual PetscErrorCode solvePoisson();

    /** \brief Project the velocity field onto the divergence-free space and
     * update the pressure field. */
    virtual PetscErrorCode projectionStep();

    /** \brief Create operators. */
    virtual PetscErrorCode createOperators();

    /** \brief Create vectors. */
    virtual PetscErrorCode createVectors();

    /** \brief Set Poisson nullspace or pin pressure at a reference point. */
    virtual PetscErrorCode setNullSpace();

    /** \brief Create an ASCII PetscViewer.
     *
     * \param filePath [in] Path of the file to write in
     * \param mode [in] File mode
     * \param viewer [out] PetscViewer object
     * \return PetscErrorCode
     */
    PetscErrorCode createPetscViewerASCII(const std::string &filePath,
                                          const PetscFileMode &mode,
                                          PetscViewer &viewer);

    /** \brief Write the solution fields into a HDF5 file.
     *
     * \param filePath [in] Path of the file to write in
     * \return PetscErrorCode
     */
    PetscErrorCode writeSolutionHDF5(const std::string &filePath);

    /** \brief Write data required to restart a simulation into a HDF5 file.
     *
     * \param filePath [in] Path of the file to write in
     * \return PetscErrorCode
     */
    virtual PetscErrorCode writeRestartDataHDF5(const std::string &filePath);

    /** \brief Read data required to restart a simulation from a HDF5 file.
     *
     * \param filePath [in] Path of the file to read from
     * \return PetscErrorCode
     */
    virtual PetscErrorCode readRestartDataHDF5(const std::string &filePath);

    /** \brief Write numbers of iterations and residuals of solvers to file. */
    virtual PetscErrorCode writeLinSolversInfo();

    /** \brief Write the time value into a HDF5 file.
     *
     * \param t [in] Time
     * \param filePath [in] Path of the file to write in
     * \return PetscErrorCode
     */
    PetscErrorCode writeTimeHDF5(const PetscReal &t,
                                 const std::string &filePath);

    /** \brief Read the time value from a HDF5 file.
     *
     * \param filePath [in] Path of the file to read from
     * \param t [out] Time
     * \return PetscErrorCode
     */
    PetscErrorCode readTimeHDF5(const std::string &filePath,
                                PetscReal &t);

    /** \brief Monitor the solution at probes. */
    PetscErrorCode monitorProbes();

};  // NavierStokesSolver
