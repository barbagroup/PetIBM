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

#include <petibm/type.h>
#include <petibm/mesh.h>
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

    /** \brief MPI communicator. */
    MPI_Comm comm;

    /** \brief Number of processes in the MPI communicator. */
    PetscMPIInt commSize;

    /** \brief Rank of the local process in the MPI communicator. */
    PetscMPIInt commRank;

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

    /** \brief Setup the probe.
     *
     * Create a sub-mesh, write the sub-mesh coordinates,
     * and create the index set.
     *
     * \param mesh [in] Cartesian mesh object
     * \return PetscErrorCode
     */
    virtual PetscErrorCode setUp(const type::Mesh &mesh);

    /** \brief Create the index set for the points to monitor.
     *
     * \param mesh [in] Cartesian mesh object
     * \return PetscErrorCode
     */
    virtual PetscErrorCode createIS(const type::Mesh &mesh);

    /** \brief Get the sub-mesh area to monitor.
     * 
     * \param mesh [in] Cartesian mesh object
     * \return PetscErrorCode
     */
    virtual PetscErrorCode createSubMesh(const type::Mesh &mesh);

    /** \brief Write the sub mesh grid points into a file.
     *
     * \param filePath [in] Path of the file to write in
     * \return PetscErrorCode
     */
    virtual PetscErrorCode writeSubMesh(const std::string &filePath);
    
    /** \brief Write the sub vector data into a file.
     *
     * \param fvec [in] Full-vector of the field solution
     * \param t [in] Time
     * \return PetscErrorCode
     */
    virtual PetscErrorCode writeSubVec(const Vec &fvec, const PetscReal &t);

    /** \brief Monitor the field solution at points in the index set.
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

    /** \brief Initialize the probe.
     * 
     * \param comm [in] MPI communicator
     * \param node [in] YAML configuration node
     * \return PetscErrorCode
     */
    virtual PetscErrorCode init(const MPI_Comm &comm,
                                const YAML::Node &node,
                                const type::Mesh &mesh) = 0;

    /** \brief Get information about the sub-mesh area to monitor.
     *
     * \param mesh [in] Cartesian mesh object
     * \param box [in] Box area to monitor
     * \return PetscErrorCode
     */
    PetscErrorCode getSubMeshInfo(const type::Mesh &mesh,
                                  const type::RealVec2D &box);

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

    /** \brief Write the sub-vector data into an ASCII file.
     * 
     * \param fvec [in] Full-vector with field data
     * \param t [in] Time
     * \return PetscErrorCode
     */
    PetscErrorCode writeSubVecASCII(const Vec &fvec, const PetscReal &t);

    /** \brief Write the sub-vector data into a HDF5 file.
     * 
     * \param fvec [in] Full-vector with field data
     * \param t [in] Time
     * \return PetscErrorCode
     */
    PetscErrorCode writeSubVecHDF5(const Vec &fvec, const PetscReal &t);

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

protected:
    /** \brief Limits of the volume. */
    type::RealVec2D box;

    /** \copydoc ProbeBase::init() */
    PetscErrorCode init(const MPI_Comm &comm,
                        const YAML::Node &node,
                        const type::Mesh &mesh);

};  // ProbeVolume

/**
 * \brief Probe class to monitor a minimal volume around a point.
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

protected:
    /** \brief Coordinates of the point to monitor around. */
    type::RealVec1D loc;

    /** \copydoc ProbeBase::init() */
    PetscErrorCode init(const MPI_Comm &comm,
                        const YAML::Node &node,
                        const type::Mesh &mesh);

    /** \brief Get the box surrounding a given point.
     *
     * \param mesh [in] Cartesian mesh object
     * \param loc [in] Coordinates of the point
     * \param box [in] Box surrounding the point
     * \return PetscErrorCode
     */
    PetscErrorCode getBox(const type::Mesh &mesh,
                          const type::RealVec1D &loc,
                          type::RealVec2D &box);

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
