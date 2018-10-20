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

#include <petibm/interp.h>
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
    /** \brief Name of the probe as a string. */
    std::string name;

    /** \brief Path of the file to output the solution. */
    std::string path;

    /** \brief Type of the field to monitor. */
    type::Field field;

    /** \brief Frequency of saving (number of time steps). */
    PetscInt nsave;

    /** \brief Monitoring starting time. */
    PetscReal tstart;

    /** \brief Monitoring ending time. */
    PetscReal tend;

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
     * \param t [in] Time
     * \return PetscErrorCode
     */
    virtual PetscErrorCode monitor(const type::Solution &solution,
                                   const type::Mesh &mesh,
                                   const PetscReal &t);

protected:
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
    
    virtual PetscErrorCode writeData(const type::Mesh &mesh,
                                     const Vec &fvec,
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
    Vec vec;

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

    PetscErrorCode writeData(const type::Mesh &mesh,
                             const Vec &fvec,
                             const PetscReal &t);

    /** \brief Write the sub vector data into a file.
     *
     * \param t [in] Time
     * \return PetscErrorCode
     */
    PetscErrorCode writeSubVec(const PetscReal &t);

    /** \brief Write the sub-vector data into an ASCII file.
     * 
     * \param t [in] Time
     * \return PetscErrorCode
     */
    PetscErrorCode writeSubVecASCII(const PetscReal &t);

    /** \brief Write the sub-vector data into a HDF5 file.
     * 
     * \param t [in] Time
     * \return PetscErrorCode
     */
    PetscErrorCode writeSubVecHDF5(const PetscReal &t);

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

    /** \brief Interpolated value. */
    PetscReal value;

    /** \brief Interpolating object to monitor at a single point. */
    type::LinInterp interp;

    /** \copydoc ProbeBase::init() */
    PetscErrorCode init(const MPI_Comm &comm,
                        const YAML::Node &node,
                        const type::Mesh &mesh);

    PetscErrorCode writeData(const type::Mesh &mesh,
                             const Vec &fvec,
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
