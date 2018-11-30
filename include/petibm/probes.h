/**
 * \file probes.h
 * \brief Prototype of probe classes, type::Probe, and factory function.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 */

#pragma once

#include <petscis.h>
#include <petscsys.h>

#include <yaml-cpp/yaml.h>

#include <petibm/lininterp.h>
#include <petibm/mesh.h>
#include <petibm/type.h>
#include <petibm/solution.h>

namespace petibm
{
namespace misc
{
/**
 * \brief Abstract Base Class of a probe.
 *
 * \see miscModule, petibm::type::Probe, petibm::misc::createProbe
 * \ingroup miscModule
 */
class ProbeBase
{
public:
    /** \brief Default constructor. */
    ProbeBase() = default;

    /** \brief Constructor. Initialize the probe.
     *
     * \param comm [in] MPI communicator
     * \param node [in] YAML configuration node
     */
    ProbeBase(const MPI_Comm &comm,
              const YAML::Node &node,
              const type::Mesh &mesh);

    /** \brief Destructor. */
    virtual ~ProbeBase();

    /** \brief Manually destroy data in the object. */
    virtual PetscErrorCode destroy();

    /** \brief Monitor the field solution and output data to file.
     *
     * \param solution [in] Data object with the field solutions
     * \param mesh [in] Cartesian mesh object
     * \param n [in] Time-step index
     * \param t [in] Time
     * \return PetscErrorCode
     */
    PetscErrorCode monitor(const type::Solution &solution,
                           const type::Mesh &mesh,
                           const PetscInt &n,
                           const PetscReal &t);

protected:
    /** \brief Name of the probe as a string. */
    std::string name;

    /** \brief Path of the file to output the solution. */
    std::string path;

    /** \brief Type of the field to monitor. */
    type::Field field;

    /** \brief Frequency of monitoring the solution (number of time steps). */
    PetscInt n_monitor;

    /** \brief Monitoring starting time. */
    PetscReal t_start;

    /** \brief Monitoring ending time. */
    PetscReal t_end;

    /** \brief PETSc viewer to output the solution. */
    PetscViewer viewer;

    /** \brief Type of the PETSc viewer to use. */
    PetscViewerType viewerType;

    /** \brief MPI communicator. */
    MPI_Comm comm;

    /** \brief Number of processes in the MPI communicator. */
    PetscMPIInt commSize;

    /** \brief Rank of the local process in the MPI communicator. */
    PetscMPIInt commRank;

    /** \brief Initialize the probe.
     * 
     * \param comm [in] MPI communicator
     * \param node [in] YAML configuration node
     * \return PetscErrorCode
     */
    virtual PetscErrorCode init(const MPI_Comm &comm,
                                const YAML::Node &node,
                                const type::Mesh &mesh) = 0;
    
    /** \brief Monitor a sub-region of a vector
     *         and possibly output data to file.
     *
     * \param da [in] Parallel layout of the full-domain vector
     * \param fvec [in] Full-domain PETSc Vec object to monitor
     * \param n [in] Time-step index
     * \param t [in] Time
     */
    virtual PetscErrorCode monitorVec(const DM &da,
                                      const Vec &fvec,
                                      const PetscInt &n,
                                      const PetscReal &t) = 0;

};  // ProbeBase

/**
 * \brief Probe class to monitor a volume region of the domain.
 *
 * \see miscModule, petibm::type::Probe, petibm::misc::createProbe
 * \ingroup miscModule
 */
class ProbeVolume : public ProbeBase
{
public:
    /** \copydoc ProbeBase::ProbeBase() */
    ProbeVolume(const MPI_Comm &comm,
                const YAML::Node &node,
                const type::Mesh &mesh);

    /** \brief Default destructor. */
    virtual ~ProbeVolume() = default;

    /** \brief Manually destroy the data. */
    PetscErrorCode destroy();

protected:
    /** \brief Limits of the volume. */
    type::RealVec2D box;

    /** \brief Index set for the grid points to monitor. */
    IS is;

    /** \brief Sub-vector of the region to monitor. */
    Vec dvec;

    /** \brief Grid point coordinates in the volume. */
    type::RealVec2D coord;

    /** \brief Number of grid points along each direction in the volume. */
    type::IntVec1D nPtsDir;

    /** \brief Number of grid points in the volume. */
    PetscInt nPts;

    /** \brief Index of the first point in the volume in each direction. */
    type::IntVec1D startIdxDir;

    /** \brief Absolute tolerance criterion when comparing values. */
    PetscReal atol;

    /** \brief Frequency of saving the data to file. */
    PetscInt n_sum;

    /** \brief Counter to know when to flush to the data to file. */
    PetscInt count = 0;

    /** \copydoc ProbeBase::init() */
    PetscErrorCode init(const MPI_Comm &comm,
                        const YAML::Node &node,
                        const type::Mesh &mesh);

    /** \brief Get information about the sub-mesh area to monitor.
     *
     * \param mesh [in] Cartesian mesh object
     * \param box [in] Box area to monitor
     * \return PetscErrorCode
     */
    PetscErrorCode getSubMeshInfo(const type::Mesh &mesh,
                                  const type::RealVec2D &box);

    /** \brief Create the index set for the points to monitor.
     *
     * \param mesh [in] Cartesian mesh object
     * \return PetscErrorCode
     */
    PetscErrorCode createIS(const type::Mesh &mesh);

    /** \brief Get the sub-mesh area to monitor.
     * 
     * \param mesh [in] Cartesian mesh object
     * \return PetscErrorCode
     */
    PetscErrorCode createSubMesh(const type::Mesh &mesh);

    /** \brief Write the sub mesh grid points into a file.
     *
     * \param filePath [in] Path of the file to write in
     * \return PetscErrorCode
     */
    PetscErrorCode writeSubMesh(const std::string &filePath);

    /** \brief Write the sub mesh into an ASCII file.
     *
     * \param filePath [in] Path of the file to write in
     * \return PetscErrorCode
     */
    PetscErrorCode writeSubMeshASCII(const std::string &filePath);

    /** \brief Write the sub mesh into a HDF5 file.
     *
     * \param filePath [in] Path of the file to write in
     * \return PetscErrorCode
     */
    PetscErrorCode writeSubMeshHDF5(const std::string &filePath);

    /** \copydoc ProbeBase::monitorVec() */
    PetscErrorCode monitorVec(const DM &da,
                              const Vec &fvec,
                              const PetscInt &n,
                              const PetscReal &t);

    /** \brief Output a PETSc Vec object to file.
     *
     * Supported formats are HDF5 and ASCII.
     *
     * \param vec [in] PETSc Vec object to output
     * \param t [in] Time
     * \return PetscErrorCode
     */
    PetscErrorCode writeVec(const Vec &vec, const PetscReal &t);

    /** \brief Output a PETSc Vec object to an ASCII file.
     *
     * \param vec [in] PETSc Vec object to output
     * \param t [in] Time
     * \return PetscErrorCode
     */
    PetscErrorCode writeVecASCII(const Vec &vec, const PetscReal &t);

    /** \brief Output a PETSc Vec object to a HDF5 file.
     *
     * \param vec [in] PETSc Vec object to output
     * \param t [in] Time
     * \return PetscErrorCode
     */
    PetscErrorCode writeVecHDF5(const Vec &vec, const PetscReal &t);

};  // ProbeVolume

/**
 * \brief Probe class to monitor the solution at a single point.
 *
 * \see miscModule, petibm::type::Probe, petibm::misc::createProbe
 * \ingroup miscModule
 */
class ProbePoint : public ProbeBase
{
public:
    /** \copydoc ProbeBase::ProbeBase() */
    ProbePoint(const MPI_Comm &comm,
               const YAML::Node &node,
               const type::Mesh &mesh);

    /** \brief Default destructor. */
    virtual ~ProbePoint() = default;

    /** \brief Manually destroy the data. */
    PetscErrorCode destroy();

protected:
    /** \brief Coordinates of the point to monitor around. */
    type::RealVec1D loc;

    /** \brief True if target point located on local sub-domain. */
    PetscBool pointOnLocalProc;

    /** \brief Local (ghosted) PETSc Vec object with neighboring values. */
    Vec svec;

    /** \brief Interpolated value. */
    PetscReal value;

    /** \brief Interpolating object to monitor at a single point. */
    type::LinInterp interp;

    /** \copydoc ProbeBase::init() */
    PetscErrorCode init(const MPI_Comm &comm,
                        const YAML::Node &node,
                        const type::Mesh &mesh);

    /** \copydoc ProbeBase::monitorVec() */
    PetscErrorCode monitorVec(const DM &da,
                              const Vec &fvec,
                              const PetscInt &n,
                              const PetscReal &t);

};  // ProbePoint

}  // end of namespace misc

namespace type
{
/**
 * \brief Type definition of Probe.
 *
 * Please use petibm::misc::createProbe to create a Probe object.
 *
 * Example usage:
 * \code
 * PetscErrorCode ierr;
 * petibm::misc::Probe probe;
 * petibm::mesh::Mesh mesh;
 * YAML::node config;
 *
 * ierr = petibm::misc::createProbe(PETSC_COMM_WORLD, config,
 *                                  mesh, probe); CHKERRQ(ierr);
 * \endcode
 *
 * \see miscModule, petibm::misc::createProbe, petibm::misc::ProbeBase
 * \ingroup miscModule
 */
typedef std::shared_ptr<misc::ProbeBase> Probe;

}  // end of namespace type

namespace misc
{
/**
 * \brief Factory function to create a probe to monitor the solution.
 *
 * \param comm [in] MPI communicator
 * \param node [in] YAML configuration node
 * \param mesh [in] Cartesian mesh object
 * \param probe [out] Probe
 * \return PetscErrorCode
 *
 * \see miscModule, petibm::type::Probe
 * \ingroup miscModule
 */
PetscErrorCode createProbe(const MPI_Comm &comm,
                           const YAML::Node &node,
                           const type::Mesh &mesh,
                           type::Probe &probe);

}  // end of namespace misc

}  // end of namespace petibm
