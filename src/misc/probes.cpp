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
    type::ProbeType probeType;

    PetscFunctionBeginUser;

    probeType = type::str2ProbeType[node["type"].as<std::string>("VOLUME")];

    // either create a probe to monitor a sub-volume or to monitor a single point
    switch (probeType)
    {
        case type::ProbeType::VOLUME:
            probe = std::make_shared<ProbeVolume>(comm, node, mesh);
            break;
        case type::ProbeType::POINT:
            probe = std::make_shared<ProbePoint>(comm, node, mesh);
            break;
        default:
            SETERRQ(comm, PETSC_ERR_ARG_UNKNOWN_TYPE,
                    "Unknown type of probe. Accepted values are: "
                    "VOLUME and POINT");
    }

    PetscFunctionReturn(0);
}  // createProbe

//***************************************************************************//
//*************************        ProbeBase        *************************//
//***************************************************************************//

// Initialize the probe.
PetscErrorCode ProbeBase::init(const MPI_Comm &inComm,
                               const YAML::Node &node,
                               const type::Mesh &mesh)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    // store information about the communicator
    comm = inComm;
    ierr = MPI_Comm_size(comm, &commSize); CHKERRQ(ierr);
    ierr = MPI_Comm_rank(comm, &commRank); CHKERRQ(ierr);

    // store information about what to monitor and when to monitor
    name = node["name"].as<std::string>("unnamed");
    field = type::str2fd[node["field"].as<std::string>()];
    path = node["path"].as<std::string>();
    n_monitor = node["n_monitor"].as<PetscInt>(1);
    t_start = node["t_start"].as<PetscReal>(0.0);
    t_end = node["t_end"].as<PetscReal>(1e12);

    viewer = PETSC_NULL;

    PetscFunctionReturn(0);
}  // ProbeBase::init

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
    if (viewer != PETSC_NULL) {ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);}
    comm = MPI_COMM_NULL;
    commSize = commRank = 0;

    PetscFunctionReturn(0);
}  // ProbeBase::destroy

// Monitor the field solution at points in the index set.
PetscErrorCode ProbeBase::monitor(const type::Solution &solution,
                                  const type::Mesh &mesh,
                                  const PetscInt &n,
                                  const PetscReal &t)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    // monitor only if appropriate time-step index and appropriate time
    if (n % n_monitor == 0 && t >= t_start && t <= t_end)
    {
        if (field < mesh->dim)  // monitor a component of the velocity field
        {
            std::vector<Vec> vel(mesh->dim);
            ierr = DMCompositeGetAccessArray(mesh->UPack, solution->UGlobal,
                                            mesh->dim, nullptr,
                                            vel.data()); CHKERRQ(ierr);
            ierr = monitorVec(mesh->da[field], vel[field], n, t); CHKERRQ(ierr);
            ierr = DMCompositeRestoreAccessArray(mesh->UPack, solution->UGlobal,
                                                mesh->dim, nullptr,
                                                vel.data()); CHKERRQ(ierr);
        }
        else  // monitor the pressure field
        {
            ierr = monitorVec(mesh->da[3], solution->pGlobal , n, t); CHKERRQ(ierr);
        }
    }

    PetscFunctionReturn(0);
}  // ProbeBase::monitor

//***************************************************************************//
//*************************       ProbeVolume       *************************//
//***************************************************************************//

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

    // store information about the type of PETSc Viewer object to use
    std::string vtype_str = node["viewer"].as<std::string>("ascii");
    if (vtype_str == "ascii")
        viewerType = PETSCVIEWERASCII;
    else if (vtype_str == "hdf5")
        viewerType = PETSCVIEWERHDF5;
    // tolerance to define if a point belong to the sub-volume
    atol = node["atol"].as<PetscReal>(1e-6);
    // number of the time-steps over which we accumulate the data
    // data are added together and we write the time averaged data
    n_sum = node["n_sum"].as<PetscInt>(0);

    is = PETSC_NULL;
    dvec = PETSC_NULL;

    // store information about the sub-volume to monitor
    box = type::RealVec2D(3, type::RealVec1D(2, 0.0));
    for (auto item : node["box"])
    {
        type::Dir dir = type::str2dir[item.first.as<std::string>()];
        box[dir][0] = item.second[0].as<PetscReal>();
        box[dir][1] = item.second[1].as<PetscReal>();
    }

    nPtsDir.resize(3, 1);
    startIdxDir.resize(3, 0);

    // get information about the location of the sub-mesh
    ierr = getSubMeshInfo(mesh, box); CHKERRQ(ierr);
    // create gridline coordinates for the sub-mesh
    ierr = createSubMesh(mesh); CHKERRQ(ierr);
    // write the sub-mesh to file
    ierr = writeSubMesh(path); CHKERRQ(ierr);
    // create a PETSc Index Set object to easily grab a sub-vector
    ierr = createIS(mesh); CHKERRQ(ierr);

    // create a PETSc Viewer to output the data
    ierr = PetscViewerCreate(comm, &viewer); CHKERRQ(ierr);
    ierr = PetscViewerSetType(viewer, viewerType); CHKERRQ(ierr);
    // Note: we set the "append" mode as the output file already exists
    // (it was created when the sub-mesh was written into it)
    ierr = PetscViewerFileSetMode(viewer, FILE_MODE_APPEND); CHKERRQ(ierr);
    ierr = PetscViewerFileSetName(viewer, path.c_str()); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // ProbeVolume::init

// Manually destroy data in the object.
PetscErrorCode ProbeVolume::destroy()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = ProbeBase::destroy(); CHKERRQ(ierr);
    if (is != PETSC_NULL) {ierr = ISDestroy(&is); CHKERRQ(ierr);}
    if (dvec != PETSC_NULL) {ierr = VecDestroy(&dvec); CHKERRQ(ierr);}

    PetscFunctionReturn(0);
}  // ProbeVolume::destroy

// Get information about the sub-mesh area to monitor.
PetscErrorCode ProbeVolume::getSubMeshInfo(const type::Mesh &mesh,
                                           const type::RealVec2D &box)
{
    std::vector<PetscReal>::iterator low, up;

    PetscFunctionBeginUser;

    // get the starting index along a gridline and the number of points
    // for each direction of the sub-mesh
    for (PetscInt d = 0; d < mesh->dim; ++d)
    {
        type::RealVec1D line(mesh->coord[field][d],
                             mesh->coord[field][d] + mesh->n[field][d]);
        low = std::lower_bound(line.begin(), line.end(), box[d][0] - atol);
        up = std::upper_bound(line.begin(), line.end(), box[d][1] + atol);
        startIdxDir[d] = low - line.begin();
        nPtsDir[d] = up - line.begin() - startIdxDir[d];
    }

    // get the number of points in the sub-volume
    nPts = std::accumulate(nPtsDir.begin(), nPtsDir.end(),
                           1, std::multiplies<PetscInt>());

    PetscFunctionReturn(0);
}  // ProbeVolume::getSubMeshInfo

// Create the index set for the points to monitor.
PetscErrorCode ProbeVolume::createIS(const type::Mesh &mesh)
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
}  // ProbeVolume::createIS

// Get the sub-mesh area to monitor.
PetscErrorCode ProbeVolume::createSubMesh(const type::Mesh &mesh)
{
    PetscFunctionBeginUser;

    // store the gridline coordinates of the sub-mesh
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
}  // ProbeVolume::createSubMesh

// Write the sub-mesh grid points into a file.
PetscErrorCode ProbeVolume::writeSubMesh(const std::string &filePath)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    // only ASCII and HDF5 formats are currently supported
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
}  // ProbeVolume::writeSubMesh

// Write the sub mesh into a HDF5 file.
PetscErrorCode ProbeVolume::writeSubMeshHDF5(const std::string &filePath)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    // only the first process in the communicator write the sub-mesh into a file
    if (commRank == 0)
    {
        // because only one process is involved in writing the sub-mesh,
        // we need to create a temporary viewer
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

    PetscFunctionReturn(0);
}  // ProbeVolume::writeSubMeshHDF5

// Write the sub mesh into an ASCII file.
PetscErrorCode ProbeVolume::writeSubMeshASCII(const std::string &filePath)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    // only the first process in the communicator write the sub-mesh into a file
    if (commRank == 0)
    {
        // because only one process is involved in writing the sub-mesh,
        // we need to create a temporary viewer
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
}  // ProbeVolume::writeSubMeshASCII

// Monitor a sub-region of the full-domain PETSc Vec object.
PetscErrorCode ProbeVolume::monitorVec(const DM &da,
                                       const Vec &fvec,
                                       const PetscInt &n,
                                       const PetscReal &t)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    // grab the part of the vector that corresponds to the sub-volume
    Vec svec;
    ierr = VecGetSubVector(fvec, is, &svec); CHKERRQ(ierr);
    if (n_sum != 0)  // we accumulate the data over the time-steps
    {
        if (dvec == PETSC_NULL)
        {
            ierr = VecDuplicate(svec, &dvec); CHKERRQ(ierr);
            ierr = VecSet(dvec, 0.0); CHKERRQ(ierr);
        }
        ierr = VecAXPY(dvec, 1.0, svec); CHKERRQ(ierr);
        count++;
        if (count % n_sum == 0) // output the time-averaged data
        {
            ierr = VecScale(dvec, 1.0 / count); CHKERRQ(ierr);
            ierr = writeVec(dvec, t); CHKERRQ(ierr);
            // reset for the next accumulation cycle
            ierr = VecSet(dvec, 0.0); CHKERRQ(ierr);
            count = 0;
        }
    }
    else  // output the time-step sub-volume data to file
    {
        ierr = writeVec(svec, t); CHKERRQ(ierr);
    }
    ierr = VecRestoreSubVector(fvec, is, &svec); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // ProbeVolume::monitorVec

// Write vector data into a file.
PetscErrorCode ProbeVolume::writeVec(const Vec &vec, const PetscReal &t)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    // only ASCII and HDF5 formats are currently supported
    if (std::strcmp(viewerType, "ascii") == 0)
    {
        ierr = writeVecASCII(vec, t); CHKERRQ(ierr);
    }
    else if (std::strcmp(viewerType, "hdf5") == 0)
    {
        ierr = writeVecHDF5(vec, t); CHKERRQ(ierr);
    }
    else
        SETERRQ(comm, PETSC_ERR_ARG_UNKNOWN_TYPE,
                "Unsupported PetscViewerType. Supported types are:\n"
                "\t PETSCVIEWERASCII and PETSCVIEWERHDF5");

    PetscFunctionReturn(0);
}  // ProbeVolume::writeVec

// Write vector data data into a HDF5 file.
PetscErrorCode ProbeVolume::writeVecHDF5(const Vec &vec, const PetscReal &t)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    std::string field_str = type::fd2str[field];
    ierr = PetscViewerHDF5PushGroup(viewer, field_str.c_str()); CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) vec,
                              std::to_string(t).c_str()); CHKERRQ(ierr);
    ierr = VecView(vec, viewer); CHKERRQ(ierr);
    if (count != 0)  // add number of time-steps accumulated as attribute
    {
        ierr = PetscViewerHDF5WriteAttribute(
            viewer, ("/" + field_str + "/" + std::to_string(t)).c_str(), "count",
            PETSC_INT, &count); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}  // ProbeVolume::writeVecHDF5

// Write vector data into an ASCII file.
PetscErrorCode ProbeVolume::writeVecASCII(const Vec &vec, const PetscReal &t)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_SYMMODU); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer, "\nt = %e\n", t); CHKERRQ(ierr);
    if (count != 0)  // write number of time-steps accumulated
    {
        ierr = PetscViewerASCIIPrintf(viewer, "count = %d\n", count); CHKERRQ(ierr);
    }
    ierr = VecView(vec, viewer); CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(viewer); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // ProbeVolume::writeVecASCII

//***************************************************************************//
//*************************       ProbePoint        *************************//
//***************************************************************************//

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

    // only ASCII output format is supported
    viewerType = PETSCVIEWERASCII;

    // get location of the point to interpolate
    loc = type::RealVec1D(3, 0.0);
    for (PetscInt d = 0; d < mesh->dim; ++d)
        loc[d] = node["loc"][d].as<PetscReal>();

    // is the target point in on the current process sub-domain?
    pointOnLocalProc = mesh->isPointOnLocalProc(loc, field);

    // create a PETSc Viewer to output the data
    // Note: as the target point belongs to only one sub-domain,
    // we create a Viewer that involves only one process
    ierr = PetscViewerCreate(PETSC_COMM_SELF, &viewer); CHKERRQ(ierr);
    ierr = PetscViewerSetType(viewer, viewerType); CHKERRQ(ierr);
    ierr = PetscViewerFileSetMode(viewer, FILE_MODE_WRITE); CHKERRQ(ierr);
    ierr = PetscViewerFileSetName(viewer, path.c_str()); CHKERRQ(ierr);
    if (pointOnLocalProc)
    {
        // create the interpolation object (2D: bi-linear, 3D: tri-linear)
        ierr = createLinInterp(PETSC_COMM_SELF,
                               loc, mesh, field, interp); CHKERRQ(ierr);
    }
    // create a local vector that will have ghost-point values
    ierr = DMCreateLocalVector(mesh->da[field], &svec); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // ProbePoint::init

// Manually destroy data in the object.
PetscErrorCode ProbePoint::destroy()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = ProbeBase::destroy(); CHKERRQ(ierr);
    if (pointOnLocalProc)
    {
        ierr = interp->destroy(); CHKERRQ(ierr);
    }
    ierr = VecDestroy(&svec); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // ProbePoint::destroy

// Monitor a sub-region of the full-domain PETSc Vec object.
PetscErrorCode ProbePoint::monitorVec(const DM &da,
                                      const Vec &fvec,
                                      const PetscInt &n,
                                      const PetscReal &t)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    // scatter values to local vector
    ierr = DMGlobalToLocalBegin(da, fvec, INSERT_VALUES, svec); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(da, fvec, INSERT_VALUES, svec); CHKERRQ(ierr);

    if (pointOnLocalProc)
    {
        // interpolate the value at the target point and write it into file
        ierr = interp->interpolate(da, svec, value); CHKERRQ(ierr);
        ierr = PetscViewerASCIIPrintf(viewer, "%10.8e\t%10.8e\n", t, value); CHKERRQ(ierr);
        ierr = PetscViewerFileSetMode(viewer, FILE_MODE_APPEND); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}  // ProbePoint::monitorVec

}  // end of namespace misc

}  // end of namespace petibm
