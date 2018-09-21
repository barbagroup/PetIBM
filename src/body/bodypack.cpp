/**
 * \file bodypack.cpp
 * \brief Implementations of body::BodyPackBase and factory function.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 */

#include <petibm/bodypack.h>
#include <petibm/io.h>

namespace petibm
{
namespace body
{
BodyPackBase::BodyPackBase(const MPI_Comm &comm, const PetscInt &dim,
                           const YAML::Node &node)
{
    init(comm, dim, node);
}  // BodyPackBase

BodyPackBase::~BodyPackBase()
{
    PetscErrorCode ierr;
    PetscBool finalized;

    PetscFunctionBeginUser;

    ierr = PetscFinalized(&finalized); CHKERRV(ierr);
    if (finalized) return;

    ierr = DMDestroy(&dmPack); CHKERRV(ierr);
    comm = MPI_COMM_NULL;
}  // ~BodyPackBase

PetscErrorCode BodyPackBase::destroy()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    dim = -1;
    nBodies = nPts = nLclPts = 0;
    std::vector<type::SingleBody>().swap(bodies);
    ierr = DMDestroy(&dmPack); CHKERRQ(ierr);
    info = "";

    comm = MPI_COMM_NULL;
    mpiSize = mpiRank = 0;
    type::IntVec1D().swap(nLclAllProcs);
    type::IntVec1D().swap(offsetsAllProcs);

    PetscFunctionReturn(0);
}  // destroy

PetscErrorCode BodyPackBase::init(const MPI_Comm &inComm, const PetscInt &inDim,
                                  const YAML::Node &node)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    // get info from MPI communicator
    comm = inComm;
    ierr = MPI_Comm_size(comm, &mpiSize); CHKERRQ(ierr);
    ierr = MPI_Comm_rank(comm, &mpiRank); CHKERRQ(ierr);

    // set the dimension
    dim = inDim;

    // get the number of bodies
    if (!node["bodies"].IsDefined())
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                "No key \"bodies\" found in the YAML node provided for "
                "initializing a BodyPack instance.\n");

    nBodies = node["bodies"].size();

    // re-size the vector holding all SingleBody instances
    bodies.resize(nBodies);

    // allocate other vectors
    nLclAllProcs = type::IntVec1D(mpiSize, 0);
    offsetsAllProcs = type::IntVec1D(mpiSize, 0);

    // initialize variables
    nPts = 0;
    nLclPts = 0;

    // loop through all bodies in the YAML node
    for (PetscInt i = 0; i < nBodies; ++i)
    {
        std::string name, filePath, type;

        // if user didn't set the name of body, use the index in vector `bodies`
        name = node["bodies"][i]["name"].as<std::string>("body" +
                                                         std::to_string(i));

        // so far, we only have one type: points. So this is also the default
        // value
        type = node["bodies"][i]["type"].as<std::string>("points");

        // check if the key "file" exists
        if (!node["bodies"][i]["file"].IsDefined())
            SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                     "No key \"file\" found in the YAML node of the body "
                     "\"%s\".\n",
                     name.c_str());

        filePath = node["bodies"][i]["file"].as<std::string>();

        // check if file path is absolute; if not, prepend directory path
        // note: this only works on Unix-like OS
        if (filePath[0] != '/')
            filePath = node["directory"].as<std::string>() + "/" + filePath;

        ierr = createSingleBody(comm, dim, type, name, filePath, bodies[i]);
        CHKERRQ(ierr);

        for (PetscMPIInt r = 0; r < mpiSize; ++r)
            nLclAllProcs[r] += bodies[i]->nLclAllProcs[r];

        // add to total number of points
        nPts += bodies[i]->nPts;

        // add to total number of local points
        nLclPts += bodies[i]->nLclPts;
    }

    // calculate offsets
    for (PetscMPIInt r = mpiSize - 1; r > 0; r--)
        offsetsAllProcs[r] = nLclAllProcs[r - 1];
    for (PetscMPIInt r = 1; r < mpiSize; r++)
        offsetsAllProcs[r] += offsetsAllProcs[r - 1];

    // create the DMComposite object (of 1D DMDA objects)
    ierr = createDmPack(); CHKERRQ(ierr);

    // create info
    ierr = createInfoString(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // init

PetscErrorCode BodyPackBase::createDmPack()
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = DMCompositeCreate(comm, &dmPack); CHKERRQ(ierr);

    for (PetscInt i = 0; i < nBodies; ++i)
    {
        ierr = DMCompositeAddDM(dmPack, bodies[i]->da); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}  // createDmPack

PetscErrorCode BodyPackBase::updateMeshIdx(const type::Mesh &mesh)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    for (auto body : bodies)
    {
        ierr = body->updateMeshIdx(mesh); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}  // updateMeshIdx

PetscErrorCode BodyPackBase::createInfoString()
{
    PetscFunctionBeginUser;

    if (mpiRank == 0)
    {
        std::stringstream ss;
        ss << std::string(80, '=') << std::endl;
        ss << "Body Pack:" << std::endl;
        ss << std::string(80, '=') << std::endl;
        ss << "\tDimension: " << dim << std::endl << std::endl;
        ss << "\tNumber of bodies: " << nBodies << std::endl << std::endl;
        ss << "\tName of bodies: [";
        for (PetscInt i = 0; i < nBodies; ++i) ss << bodies[i]->name << ", ";
        ss << "]" << std::endl;
        info = ss.str();
    }

    PetscFunctionReturn(0);
}  // createInfoString

PetscErrorCode BodyPackBase::printInfo() const
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = io::print(info); CHKERRQ(ierr);

    for (PetscInt i = 0; i < nBodies; ++i)
    {
        ierr = bodies[i]->printInfo(); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}  // printInfo

PetscErrorCode BodyPackBase::findProc(const PetscInt &bIdx,
                                      const PetscInt &ptIdx,
                                      PetscMPIInt &proc) const
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    if ((bIdx < 0) || (bIdx >= nBodies))
        SETERRQ2(comm, PETSC_ERR_ARG_SIZ,
                 "Body index %d is out of range. Total number of bodies is %d.",
                 bIdx, nBodies);

    ierr = bodies[bIdx]->findProc(ptIdx, proc); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // findProc

PetscErrorCode BodyPackBase::getGlobalIndex(const PetscInt &bIdx,
                                            const PetscInt &ptIdx,
                                            const PetscInt &dof,
                                            PetscInt &idx) const
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    if ((bIdx < 0) || (bIdx >= nBodies))
        SETERRQ2(comm, PETSC_ERR_ARG_SIZ,
                 "Body index %d is out of range. Total number of bodies is %d.",
                 bIdx, nBodies);

    ierr = bodies[bIdx]->getGlobalIndex(ptIdx, dof, idx); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // getGlobalIndex

PetscErrorCode BodyPackBase::getGlobalIndex(const PetscInt &bIdx,
                                            const MatStencil &s,
                                            PetscInt &idx) const
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = getGlobalIndex(bIdx, s.i, s.c, idx); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // getGlobalIndex

PetscErrorCode BodyPackBase::getPackedGlobalIndex(const PetscInt &bIdx,
                                                  const PetscInt &ptIdx,
                                                  const PetscInt &dof,
                                                  PetscInt &idx) const
{
    PetscErrorCode ierr;
    PetscMPIInt p;
    PetscInt unPackedIdx;

    PetscFunctionBeginUser;

    ierr = findProc(bIdx, ptIdx, p); CHKERRQ(ierr);

    ierr = getGlobalIndex(bIdx, ptIdx, dof, unPackedIdx); CHKERRQ(ierr);

    idx = offsetsAllProcs[p];

    for (PetscInt b = 0; b < bIdx; ++b) idx += bodies[b]->nLclAllProcs[p];

    idx += (unPackedIdx - bodies[bIdx]->offsetsAllProcs[p]);

    PetscFunctionReturn(0);
}  // getPackedGlobalIndex

PetscErrorCode BodyPackBase::getPackedGlobalIndex(const PetscInt &bIdx,
                                                  const MatStencil &s,
                                                  PetscInt &idx) const
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = getPackedGlobalIndex(bIdx, s.i, s.c, idx); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // getPackedGlobalIndex

PetscErrorCode BodyPackBase::calculateAvgForces(const Vec &f,
                                                type::RealVec2D &fAvg)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    std::vector<Vec> unPacked(nBodies);

    fAvg.resize(nBodies);

    ierr =
        DMCompositeGetAccessArray(dmPack, f, nBodies, nullptr, unPacked.data());
    CHKERRQ(ierr);

    for (PetscInt i = 0; i < nBodies; ++i)
    {
        fAvg[i].resize(dim);
        ierr = bodies[i]->calculateAvgForces(unPacked[i], fAvg[i]);
        CHKERRQ(ierr);
    }

    ierr = DMCompositeRestoreAccessArray(dmPack, f, nBodies, nullptr,
                                         unPacked.data()); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // calculateAvgForces

PetscErrorCode createBodyPack(const MPI_Comm &comm, const PetscInt &dim,
                              const YAML::Node &node, type::BodyPack &bodies)
{
    PetscFunctionBeginUser;

    bodies = std::make_shared<BodyPackBase>(comm, dim, node);

    PetscFunctionReturn(0);
}  // createBodyPack

}  // end of namespace body

}  // end of namespace petibm
