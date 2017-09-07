/***************************************************************************//**
 * \file BodyPack.cpp
 * \author Anush Krishnan (anus@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Definition of the members of class `BodyPack`.
 */


// PetIBM
# include "BodyPack.h"


namespace petibm
{
namespace utilities
{

/** \copydoc BodyPack::BodyPack(). */
BodyPack::BodyPack() = default;


/** \copydoc BodyPack::~BodyPack(). */
BodyPack::~BodyPack() = default;


/** \copydoc BodyPack::BodyPack(const CartesianMesh &, const YAML::Node &). */
BodyPack::BodyPack(const CartesianMesh &_mesh, const YAML::Node &node)
{
    init(_mesh, node);
}


/** \copydoc BodyPack::init. */
PetscErrorCode BodyPack::init(
        const CartesianMesh &_mesh, const YAML::Node &node)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    // save the address of input mesh instance
    // note: a bad practice for shared_ptr
    mesh = std::shared_ptr<const CartesianMesh>(&_mesh, [](const CartesianMesh*){});

    // save MPI information from the mesh
    comm = mesh->comm;
    mpiSize = mesh->mpiSize;
    mpiRank = mesh->mpiRank;

    // set the dimension
    dim = mesh->dim;

    // get the number of bodies
    nBodies = node.size();

    // sizing the vector holding all SingleBody instances
    bodies.resize(nBodies);

    // allocate other vectors
    nLclAllProcs = types::IntVec1D(mpiSize, 0);
    offsetsAllProcs = types::IntVec1D(mpiSize, 0);

    // initialize variables
    nPts = 0;
    nLclPts = 0;

    // loop through all bodies in the YAML node
    for(PetscInt i=0; i<nBodies; ++i)
    {
        std::string     name, file, type;

        // if user didn't set the name of body, use the index in vector `bodies`
        name = node[i]["name"].as<std::string>("body" + std::to_string(i));

        // TODO: should we check if user really set the key "file"?
        file = node[i]["pointsFile"].as<std::string>();

        // so far, we only have one type: points. So this is also the default value
        type = node[i]["type"].as<std::string>("points");

        ierr = bodies[i].init(*mesh, file, name); CHKERRQ(ierr);

        for(PetscMPIInt r=0; r<mpiSize; ++r)
            nLclAllProcs[r] += bodies[i].nLclAllProcs[r];

        // add to total number of points
        nPts += bodies[i].nPts;

        // add to total number of local points
        nLclPts += bodies[i].nLclPts;
    }

    // calculate offsets
    for(PetscMPIInt r=mpiSize-1; r>0; r--) offsetsAllProcs[r] = nLclAllProcs[r-1];
    for(PetscMPIInt r=1; r<mpiSize; r++) offsetsAllProcs[r] += offsetsAllProcs[r-1];

    // create dmPack
    ierr = createDmPack(); CHKERRQ(ierr);

    // create mapping
    ierr = createMappings(); CHKERRQ(ierr);

    // create info
    ierr = createInfoString(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/** \copydoc BodyPack::createDmPack. */
PetscErrorCode BodyPack::createDmPack()
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    ierr = DMCompositeCreate(*comm, &dmPack); CHKERRQ(ierr);

    for(PetscInt i=0; i<nBodies; ++i)
    {
        ierr = DMCompositeAddDM(dmPack, bodies[i].da); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}


/** \copydoc BodyPack::createMappings. */
PetscErrorCode BodyPack::createMappings()
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    ISLocalToGlobalMapping      *temp;

    ierr = DMCompositeGetISLocalToGlobalMappings(dmPack, &temp); CHKERRQ(ierr);

    mapping.assign(temp, temp+nBodies); CHKERRQ(ierr);

    // DO NOT DE-ALLOCATE THE MEMORY SPACE POINTED BY `temp`!!

    PetscFunctionReturn(0);
}


/** \copydoc BodyPack::createInfoString. */
PetscErrorCode BodyPack::createInfoString()
{
    PetscFunctionBeginUser;

    if (mpiRank == 0)
    {
        std::stringstream       ss;

        ss << std::string(80, '=') << std::endl;
        ss << "Body Pack:" << std::endl;
        ss << std::string(80, '=') << std::endl;

        ss << "\tDimension: " << dim << std::endl << std::endl;

        ss << "\tNumber of bodies: " << nBodies << std::endl << std::endl;

        ss << "\tName of bodies: " << std::endl;

        for(PetscInt i=0; i<nBodies; ++i)
            ss << "\t\t" << i << ": " << bodies[i].name << std::endl;

        ss << std::endl;

        info = ss.str();
    }

    PetscFunctionReturn(0);
}


/** \copydoc BodyPack::printInfo() const. */
PetscErrorCode BodyPack::printInfo() const
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    ierr = PetscSynchronizedPrintf(*comm, info.c_str()); CHKERRQ(ierr);
    ierr = PetscSynchronizedFlush(*comm, PETSC_STDOUT); CHKERRQ(ierr);

    for(PetscInt i=0; i<nBodies; ++i)
    {
        ierr = bodies[i].printInfo(); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}


/** \copydoc BodyPack::findProc(
 *           const PetscInt &, const PetscInt &, PetscMPIInt &) const. */
PetscErrorCode BodyPack::findProc(
        const PetscInt &bIdx, const PetscInt &ptIdx, PetscMPIInt &proc) const
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    if ((bIdx < 0) || (bIdx >= nBodies))
        SETERRQ2(*comm, PETSC_ERR_ARG_SIZ,
                "Body index %d is out of range. Total number of bodies is %d.",
                bIdx, nBodies);

    ierr = bodies[bIdx].findProc(ptIdx, proc); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/** \copydoc BodyPack::getGlobalIndex(const PetscInt &bIdx, 
        const PetscInt &ptIdx, const PetscInt &dof, PetscInt &idx) const. */
PetscErrorCode BodyPack::getGlobalIndex(const PetscInt &bIdx, 
        const PetscInt &ptIdx, const PetscInt &dof, PetscInt &idx) const
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    if ((bIdx < 0) || (bIdx >= nBodies))
        SETERRQ2(*comm, PETSC_ERR_ARG_SIZ,
                "Body index %d is out of range. Total number of bodies is %d.",
                bIdx, nBodies);

    ierr = bodies[bIdx].getGlobalIndex(ptIdx, dof, idx); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/** \copydoc BodyPack::getGlobalIndex(
        const PetscInt &bIdx, const MatStencil &s, PetscInt &idx) const. */
PetscErrorCode BodyPack::getGlobalIndex(
        const PetscInt &bIdx, const MatStencil &s, PetscInt &idx) const
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    ierr = getGlobalIndex(bIdx, s.i, s.c, idx); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/** \copydoc BodyPack::getPackedGlobalIndex(const PetscInt &bIdx, 
        const PetscInt &ptIdx, const PetscInt &dof, PetscInt &idx) const. */
PetscErrorCode BodyPack::getPackedGlobalIndex(const PetscInt &bIdx, 
        const PetscInt &ptIdx, const PetscInt &dof, PetscInt &idx) const
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    PetscMPIInt         p;
    PetscInt            unPackedIdx;

    ierr = findProc(bIdx, ptIdx, p); CHKERRQ(ierr);

    ierr = getGlobalIndex(bIdx, ptIdx, dof, unPackedIdx); CHKERRQ(ierr);

    idx = offsetsAllProcs[p];

    for(PetscInt b=0; b<bIdx; ++b)
        idx += bodies[b].nLclAllProcs[p];

    idx += (unPackedIdx - bodies[bIdx].offsetsAllProcs[p]);

    PetscFunctionReturn(0);
}


/** \copydoc BodyPack::getPackedGlobalIndex(
        const PetscInt &bIdx, const MatStencil &s, PetscInt &idx) const. */
PetscErrorCode BodyPack::getPackedGlobalIndex(
        const PetscInt &bIdx, const MatStencil &s, PetscInt &idx) const
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    ierr = getPackedGlobalIndex(bIdx, s.i, s.c, idx); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/** \copydoc BodyPack::calculateAvgForces(const Vec &, types::RealVec2D &). */
PetscErrorCode BodyPack::calculateAvgForces(
        const Vec &f, types::RealVec2D &fAvg)
{
    PetscFunctionBeginUser;

    PetscErrorCode          ierr;

    std::vector<Vec>        unPacked(nBodies);

    fAvg.resize(nBodies);

    ierr = DMCompositeGetAccessArray(
            dmPack, f, nBodies, nullptr, unPacked.data()); CHKERRQ(ierr);

    for(PetscInt i=0; i<nBodies; ++i)
    {
        fAvg[i].resize(dim);
        ierr = bodies[i].calculateAvgForces(unPacked[i], fAvg[i]); CHKERRQ(ierr);
    }

    ierr = DMCompositeRestoreAccessArray(
            dmPack, f, nBodies, nullptr, unPacked.data()); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/** \copydoc BodyPack::writeAvgForce. */
PetscErrorCode BodyPack::writeAvgForce(const PetscReal &dt,
        const Vec &f, const std::string &dir, const std::string &file)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    types::RealVec2D    fAvgs;

    PetscViewer         viewer;

    ierr = calculateAvgForces(f, fAvgs); CHKERRQ(ierr);

    ierr = PetscViewerCreate(*comm, &viewer); CHKERRQ(ierr);
    ierr = PetscViewerSetType(viewer, PETSCVIEWERASCII); CHKERRQ(ierr);
    ierr = PetscViewerFileSetMode(viewer, FILE_MODE_APPEND); CHKERRQ(ierr);
    ierr = PetscViewerFileSetName(viewer, (dir+"/"+file+".txt").c_str()); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer, "%f", dt); CHKERRQ(ierr);

    for(PetscInt i=0; i<nBodies; ++i)
    {
        for(PetscInt dof=0; dof<dim; ++dof)
        {
            ierr = PetscViewerASCIIPrintf(
                    viewer, "\t%f", fAvgs[i][dof]); CHKERRQ(ierr);
        }
    }
    ierr = PetscViewerASCIIPrintf(viewer, "\n"); CHKERRQ(ierr);

    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

} // end of namespace utilities
} // end of namespace petibm
