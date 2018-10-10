/**
 * \file probes.cpp
 * \brief Implementations of the probe classes and factory function.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 */

#include <numeric>
#include <cstring>
#include <algorithm>

#include <petscdmcomposite.h>
#include <petscviewerhdf5.h>

#include <petibm/probes.h>

namespace petibm
{

namespace misc
{
// Factory function to create a probe to monitor the solution.
PetscErrorCode createProbe(const MPI_Comm &comm,
                           const YAML::Node &node,
                           const type::Mesh &mesh,
                           type::Probe &probe)
{
    PetscErrorCode ierr;
    type::ProbeType probeType;

    PetscFunctionBeginUser;

    probeType = type::str2ProbeType[node["type"].as<std::string>("VOLUME")];

    switch (probeType)
    {
        case type::ProbeType::VOLUME:
            probe = std::make_shared<ProbeVolume>(comm, node, mesh);
            ierr = probe->setUp(mesh); CHKERRQ(ierr);
            break;
        case type::ProbeType::POINT:
            probe = std::make_shared<ProbePoint>(comm, node, mesh);
            ierr = probe->setUp(mesh); CHKERRQ(ierr);
            break;
        default:
            SETERRQ(comm, PETSC_ERR_ARG_UNKNOWN_TYPE,
                    "Unknown type of probe. Accepted values are: "
                    "VOLUME and POINT");
    }

    PetscFunctionReturn(0);
}  // createProbe

// Destructor.
ProbeBase::~ProbeBase()
{
    PetscErrorCode ierr;
    PetscBool finalized;

    PetscFunctionBeginUser;

    ierr = PetscFinalized(&finalized); CHKERRV(ierr);
    if (finalized) return;

    ierr = destroy(); CHKERRV(ierr);
}  // ProbeBase::~ProbeBase

// Manually destroy data in the object.
PetscErrorCode ProbeBase::destroy()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    name = "";
    path = "";
    if (vec != PETSC_NULL) {ierr = VecDestroy(&vec); CHKERRQ(ierr);}
    if (is != PETSC_NULL) {ierr = ISDestroy(&is); CHKERRQ(ierr);}
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
    comm = MPI_COMM_NULL;
    commSize = commRank = 0;

    PetscFunctionReturn(0);
}  // ProbeBase::destroy

// Initialize the probe.
PetscErrorCode ProbeBase::init(const MPI_Comm &inComm,
                               const YAML::Node &node,
                               const type::Mesh &mesh)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    comm = inComm;
    ierr = MPI_Comm_size(comm, &commSize); CHKERRQ(ierr);
    ierr = MPI_Comm_rank(comm, &commRank); CHKERRQ(ierr);

    name = node["name"].as<std::string>();
    field = type::str2fd[node["field"].as<std::string>()];
    path = node["path"].as<std::string>();
    nsave = node["nsave"].as<PetscInt>();
    tstart = node["tstart"].as<PetscReal>(0.0);
    tend = node["tend"].as<PetscReal>(1e12);
    std::string vtype_str = node["viewer"].as<std::string>("ascii");
    if (vtype_str == "ascii")
        viewerType = PETSCVIEWERASCII;
    else if (vtype_str == "hdf5")
        viewerType = PETSCVIEWERHDF5;
    atol = node["atol"].as<PetscReal>(1e-6);

    is = PETSC_NULL;
    vec = PETSC_NULL;

    // create a PETSc Viewer to output the data
    ierr = PetscViewerCreate(comm, &viewer); CHKERRQ(ierr);
    ierr = PetscViewerSetType(viewer, viewerType); CHKERRQ(ierr);
    ierr = PetscViewerFileSetMode(viewer, FILE_MODE_APPEND); CHKERRQ(ierr);
    ierr = PetscViewerFileSetName(viewer, path.c_str()); CHKERRQ(ierr);

    nPtsDir.resize(3, 1);
    startIdxDir.resize(3, 0);

    PetscFunctionReturn(0);
}  // ProbeBase::init

// Setup probe sub-mesh, write sub-mesh, and create index set.
PetscErrorCode ProbeBase::setUp(const type::Mesh &mesh)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = createSubMesh(mesh); CHKERRQ(ierr);
    ierr = writeSubMesh(path); CHKERRQ(ierr);
    ierr = createIS(mesh); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // setUp

// Get information about the sub-mesh area to monitor.
PetscErrorCode ProbeBase::getSubMeshInfo(const type::Mesh &mesh,
                                         const type::RealVec2D &box)
{
    std::vector<PetscReal>::iterator low, up;

    PetscFunctionBeginUser;

    for (PetscInt d = 0; d < mesh->dim; ++d)
    {
        type::RealVec1D line(mesh->coord[field][d],
                             mesh->coord[field][d] + mesh->n[field][d]);
        low = std::lower_bound(line.begin(), line.end(), box[d][0] - atol);
        up = std::upper_bound(line.begin(), line.end(), box[d][1] + atol);
        startIdxDir[d] = low - line.begin();
        nPtsDir[d] = up - line.begin() - startIdxDir[d];
    }

    // get the number of points in the volume
    nPts = std::accumulate(nPtsDir.begin(), nPtsDir.end(),
                           1, std::multiplies<PetscInt>());

    PetscFunctionReturn(0);
}  // ProbeBase::getSubMeshInfo

// Create the index set for the points to monitor.
PetscErrorCode ProbeBase::createIS(const type::Mesh &mesh)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    DMDALocalInfo info;
    ierr = DMDAGetLocalInfo(mesh->da[field], &info); CHKERRQ(ierr);
    
    std::vector<PetscInt> indices(nPts);
    PetscInt count = 0;
    for (PetscInt k = info.zs; k < info.zs + info.zm; ++k)
    {
        for (PetscInt j = info.ys; j < info.ys + info.ym; ++j)
        {
            for (PetscInt i = info.xs; i < info.xs + info.xm; ++i)
            {
                if (i >= startIdxDir[0] && i < startIdxDir[0] + nPtsDir[0] &&
                    j >= startIdxDir[1] && j < startIdxDir[1] + nPtsDir[1] &&
                    k >= startIdxDir[2] && k < startIdxDir[2] + nPtsDir[2])
                {
                    indices[count] = k * (info.my * info.mx) + j * info.mx + i;
                    count++;
                }
            }
        }
    }
    indices.resize(count);

    ierr = ISCreateGeneral(comm, count, &indices[0],
                           PETSC_COPY_VALUES, &is); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // ProbeBase::createIS

// Get the sub-mesh area to monitor.
PetscErrorCode ProbeBase::createSubMesh(const type::Mesh &mesh)
{
    PetscFunctionBeginUser;

    coord = type::RealVec2D(mesh->dim, type::RealVec1D());
    for (PetscInt d = 0; d < mesh->dim; ++d)
    {
        coord[d].reserve(nPtsDir[d]);
        PetscInt first = startIdxDir[d];
        PetscInt last = first + nPtsDir[d];
        for (PetscInt i = first; i < last; ++i)
            coord[d].push_back(mesh->coord[field][d][i]);
    }

    PetscFunctionReturn(0);
}  // ProbeBase::createSubMesh

// Monitor the field solution at points in the index set.
PetscErrorCode ProbeBase::monitor(const type::Solution &solution,
                                  const type::Mesh &mesh,
                                  const PetscReal &t)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    if (field < mesh->dim)
    {
        std::vector<Vec> vel(mesh->dim);
        ierr = DMCompositeGetAccessArray(mesh->UPack, solution->UGlobal,
                                         mesh->dim, nullptr,
                                         vel.data()); CHKERRQ(ierr);
        ierr = writeSubVec(vel[field], t); CHKERRQ(ierr);
        ierr = DMCompositeRestoreAccessArray(mesh->UPack, solution->UGlobal,
                                             mesh->dim, nullptr,
                                             vel.data()); CHKERRQ(ierr);
    }
    else
    {
        ierr = writeSubVec(solution->pGlobal, t); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}  // ProbeBase::monitor

// Write the sub mesh grid points into a file.
PetscErrorCode ProbeBase::writeSubMesh(const std::string &filePath)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    if (std::strcmp(viewerType, "ascii") == 0)
    {
        ierr = writeSubMeshASCII(filePath); CHKERRQ(ierr);
    }
    else if (std::strcmp(viewerType, "hdf5") == 0)
    {
        ierr = writeSubMeshHDF5(filePath); CHKERRQ(ierr);
    }
    else
        SETERRQ(comm, PETSC_ERR_ARG_UNKNOWN_TYPE,
                "Unsupported PetscViewerType. Supported types are:\n"
                "\t PETSCVIEWERASCII and PETSCVIEWERHDF5");

    PetscFunctionReturn(0);
}  // ProbeBase::writeSubMesh

// Write the sub vector data into a file.
PetscErrorCode ProbeBase::writeSubVec(const Vec &fvec, const PetscReal &t)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    if (std::strcmp(viewerType, "ascii") == 0)
    {
        ierr = writeSubVecASCII(fvec, t); CHKERRQ(ierr);
    }
    else if (std::strcmp(viewerType, "hdf5") == 0)
    {
        ierr = writeSubVecHDF5(fvec, t); CHKERRQ(ierr);
    }
    else
        SETERRQ(comm, PETSC_ERR_ARG_UNKNOWN_TYPE,
                "Unsupported PetscViewerType. Supported types are:\n"
                "\t PETSCVIEWERASCII and PETSCVIEWERHDF5");

    PetscFunctionReturn(0);
}  // ProbeBase::writeSubVec

// Write the sub mesh into a HDF5 file.
PetscErrorCode ProbeBase::writeSubMeshHDF5(const std::string &filePath)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    if (commRank == 0)
    {
        PetscViewer viewer2;
        ierr = PetscViewerCreate(PETSC_COMM_SELF, &viewer2); CHKERRQ(ierr);
        ierr = PetscViewerSetType(viewer2, PETSCVIEWERHDF5); CHKERRQ(ierr);
        ierr = PetscViewerFileSetMode(viewer2, FILE_MODE_WRITE); CHKERRQ(ierr);
        ierr = PetscViewerFileSetName(
            viewer2, filePath.c_str()); CHKERRQ(ierr);
        ierr = PetscViewerHDF5PushGroup(viewer2, "mesh"); CHKERRQ(ierr);
        std::vector<std::string> dirs{"x", "y", "z"};
        for (unsigned int d = 0; d < coord.size(); ++d)
        {
            Vec tmp;
            ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, 1, nPtsDir[d],
                                         &coord[d][0], &tmp); CHKERRQ(ierr);
            ierr = PetscObjectSetName((PetscObject) tmp,
                                      dirs[d].c_str()); CHKERRQ(ierr);
            ierr = VecView(tmp, viewer2); CHKERRQ(ierr);
            ierr = VecDestroy(&tmp); CHKERRQ(ierr);
        }
        ierr = PetscViewerDestroy(&viewer2); CHKERRQ(ierr);
    }
    ierr = MPI_Barrier(comm); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // ProbeBase::writeSubMeshHDF5

// Write the sub mesh into an ASCII file.
PetscErrorCode ProbeBase::writeSubMeshASCII(const std::string &filePath)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    if (commRank == 0)
    {
        PetscViewer viewer2;
        ierr = PetscViewerCreate(PETSC_COMM_SELF, &viewer2); CHKERRQ(ierr);
        ierr = PetscViewerSetType(viewer2, PETSCVIEWERASCII); CHKERRQ(ierr);
        ierr = PetscViewerPushFormat(
            viewer2, PETSC_VIEWER_ASCII_SYMMODU); CHKERRQ(ierr);
        ierr = PetscViewerFileSetMode(viewer2, FILE_MODE_WRITE); CHKERRQ(ierr);
        ierr = PetscViewerFileSetName(
            viewer2, filePath.c_str()); CHKERRQ(ierr);
        std::vector<std::string> dirs{"x", "y", "z"};
        for (unsigned int d = 0; d < coord.size(); ++d)
        {
            Vec tmp;
            ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, 1, nPtsDir[d],
                                         &coord[d][0], &tmp); CHKERRQ(ierr);
            ierr = PetscObjectSetName((PetscObject) tmp,
                                      dirs[d].c_str()); CHKERRQ(ierr);
            ierr = VecView(tmp, viewer2); CHKERRQ(ierr);
            ierr = VecDestroy(&tmp); CHKERRQ(ierr);
        }
        ierr = PetscViewerPopFormat(viewer2); CHKERRQ(ierr);
        ierr = PetscViewerDestroy(&viewer2); CHKERRQ(ierr);
    }
    ierr = MPI_Barrier(comm); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // ProbeBase::writeSubMeshASCII

// Write the sub-vector data into a HDF5 file.
PetscErrorCode ProbeBase::writeSubVecHDF5(const Vec &fvec, const PetscReal &t)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = VecGetSubVector(fvec, is, &vec); CHKERRQ(ierr);
    std::string field_str = type::fd2str[field];
    ierr = PetscViewerHDF5PushGroup(viewer, field_str.c_str()); CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) vec,
                              std::to_string(t).c_str()); CHKERRQ(ierr);
    ierr = VecView(vec, viewer); CHKERRQ(ierr);
    ierr = VecRestoreSubVector(fvec, is, &vec); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // ProbeBase::writeSubVecHDF5

// Write the sub-vector data into a ASCII file.
PetscErrorCode ProbeBase::writeSubVecASCII(const Vec &fvec, const PetscReal &t)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = VecGetSubVector(fvec, is, &vec); CHKERRQ(ierr);
    ierr = PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_SYMMODU); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer, "\nt = %e\n", t); CHKERRQ(ierr);
    ierr = VecView(vec, viewer); CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(viewer); CHKERRQ(ierr);
    ierr = VecRestoreSubVector(fvec, is, &vec); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // ProbeBase::writeSubVecASCII

// Constructor. Initialize the probe.
ProbeVolume::ProbeVolume(const MPI_Comm &comm,
                         const YAML::Node &node,
                         const type::Mesh &mesh)
{
    init(comm, node, mesh);
}  // ProbeVolume::ProbeVolume

// Initialize the probe.
PetscErrorCode ProbeVolume::init(const MPI_Comm &comm,
                                 const YAML::Node &node,
                                 const type::Mesh &mesh)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = ProbeBase::init(comm, node, mesh); CHKERRQ(ierr);

    box = type::RealVec2D(3, type::RealVec1D(2, 0.0));
    for (auto item : node["box"])
    {
        type::Dir dir = type::str2dir[item.first.as<std::string>()];
        box[dir][0] = item.second[0].as<PetscReal>();
        box[dir][1] = item.second[1].as<PetscReal>();
    }

    ierr = getSubMeshInfo(mesh, box); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // ProbeVolume::init

// Constructor. Initialize the probe.
ProbePoint::ProbePoint(const MPI_Comm &comm,
                       const YAML::Node &node,
                       const type::Mesh &mesh)
{
    init(comm, node, mesh);
}  // ProbePoint::ProbePoint

// Initialize the probe.
PetscErrorCode ProbePoint::init(const MPI_Comm &comm,
                                const YAML::Node &node,
                                const type::Mesh &mesh)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = ProbeBase::init(comm, node, mesh); CHKERRQ(ierr);

    loc = type::RealVec1D(3, 0.0);
    for (PetscInt d = 0; d < mesh->dim; ++d)
        loc[d] = node["loc"][d].as<PetscReal>();

    type::RealVec2D box = type::RealVec2D(3, type::RealVec1D(2, 0.0));
    ierr = getBox(mesh, loc, box); CHKERRQ(ierr);

    ierr = getSubMeshInfo(mesh, box); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // ProbeVolume::init

// Get the box around a given point.
PetscErrorCode ProbePoint::getBox(const type::Mesh &mesh,
                                  const type::RealVec1D &loc,
                                  type::RealVec2D &box)
{
    std::vector<PetscReal>::iterator low;

    PetscFunctionBeginUser;

    for (PetscInt d = 0; d < mesh->dim; ++d)
    {
        type::RealVec1D line(mesh->coord[field][d],
                             mesh->coord[field][d] + mesh->n[field][d]);
        low = std::lower_bound(line.begin(), line.end(), loc[d] - atol);
        PetscInt idx = low - line.begin() - 1;
        box[d][0] = line[idx] - 0.5 * (line[idx] - line[idx - 1]);
        box[d][1] = line[idx + 1] + 0.5 * (line[idx + 2] - line[idx + 1]);
    }

    PetscFunctionReturn(0);
}  // ProbePoint::getBox

}  // end of namespace misc

}  // end of namespace petibm
