/**
 * \file bodypack.cpp
 * \brief Implementations of body::BodyPackBase and factory function.
 * \author Anush Krishnan (anush@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \copyright MIT.
 */


// PetIBM
# include <petibm/bodypack.h>
# include <petibm/io.h>


namespace petibm
{
namespace body
{


BodyPackBase::BodyPackBase(const type::Mesh &inMesh, const YAML::Node &node)
{
    init(inMesh, node);
} // BodyPackBase


BodyPackBase::~BodyPackBase()
{
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    PetscBool finalized;

    ierr = PetscFinalized(&finalized); CHKERRV(ierr);
    if (finalized) return;

    ierr = DMDestroy(&dmPack); CHKERRV(ierr);
    comm = MPI_COMM_NULL;
} // ~BodyPackBase


PetscErrorCode BodyPackBase::destroy()
{
    PetscFunctionBeginUser;
    PetscErrorCode ierr;

    dim = -1;
    nBodies = nPts = nLclPts = 0;
    std::vector<type::SingleBody>().swap(bodies);
    ierr = DMDestroy(&dmPack); CHKERRQ(ierr);
    info = "";

    mesh.reset();
    comm = MPI_COMM_NULL;
    mpiSize = mpiRank = 0;
    type::IntVec1D().swap(nLclAllProcs);
    type::IntVec1D().swap(offsetsAllProcs);

    PetscFunctionReturn(0);
} // destroy


PetscErrorCode BodyPackBase::init(
        const type::Mesh &inMesh, const YAML::Node &node)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    // save the reference to background mesh
    mesh = inMesh;

    // save MPI information from the mesh
    comm = mesh->comm;
    mpiSize = mesh->mpiSize;
    mpiRank = mesh->mpiRank;

    // set the dimension
    dim = mesh->dim;

    // get the number of bodies
    if (! node["bodies"].IsDefined())
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
    for(PetscInt i=0; i<nBodies; ++i)
    {
        std::string     name, file, type;

        // if user didn't set the name of body, use the index in vector `bodies`
        name = node["bodies"][i]["name"].as<std::string>("body" + std::to_string(i));

        // so far, we only have one type: points. So this is also the default value
        type = node["bodies"][i]["type"].as<std::string>("points");

        // check if the key "file" exists
        if (! node["bodies"][i]["file"].IsDefined())
            SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                    "No key \"file\" found in the YAML node of the body "
                    "\"%s\".\n", name.c_str());
        
        file = node["bodies"][i]["file"].as<std::string>();
        
        // check if file path is absolute; if not, prepend directory path
        // note: this only works on Unix-like OS
        if (file[0] != '/') file = node["directory"].as<std::string>() + "/" + file;

        ierr = createSingleBody(mesh, type, name, file, bodies[i]); CHKERRQ(ierr);

        for(PetscMPIInt r=0; r<mpiSize; ++r)
            nLclAllProcs[r] += bodies[i]->nLclAllProcs[r];

        // add to total number of points
        nPts += bodies[i]->nPts;

        // add to total number of local points
        nLclPts += bodies[i]->nLclPts;
    }

    // calculate offsets
    for(PetscMPIInt r=mpiSize-1; r>0; r--) offsetsAllProcs[r] = nLclAllProcs[r-1];
    for(PetscMPIInt r=1; r<mpiSize; r++) offsetsAllProcs[r] += offsetsAllProcs[r-1];

    // create dmPack
    ierr = createDmPack(); CHKERRQ(ierr);

    // create info
    ierr = createInfoString(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // init


PetscErrorCode BodyPackBase::createDmPack()
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    ierr = DMCompositeCreate(comm, &dmPack); CHKERRQ(ierr);

    for(PetscInt i=0; i<nBodies; ++i)
    {
        ierr = DMCompositeAddDM(dmPack, bodies[i]->da); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
} // createDmPack


PetscErrorCode BodyPackBase::createInfoString()
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

        ss << "\tName of bodies: [";

        for(PetscInt i=0; i<nBodies; ++i) ss << bodies[i]->name << ", ";

        ss << "]" << std::endl;

        info = ss.str();
    }

    PetscFunctionReturn(0);
} // createInfoString


PetscErrorCode BodyPackBase::printInfo() const
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    ierr = io::print(info); CHKERRQ(ierr);

    for(PetscInt i=0; i<nBodies; ++i)
    {
        ierr = bodies[i]->printInfo(); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
} // printInfo


PetscErrorCode BodyPackBase::findProc(
        const PetscInt &bIdx, const PetscInt &ptIdx, PetscMPIInt &proc) const
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    if ((bIdx < 0) || (bIdx >= nBodies))
        SETERRQ2(comm, PETSC_ERR_ARG_SIZ,
                "Body index %d is out of range. Total number of bodies is %d.",
                bIdx, nBodies);

    ierr = bodies[bIdx]->findProc(ptIdx, proc); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // findProc


PetscErrorCode BodyPackBase::getGlobalIndex(const PetscInt &bIdx, 
        const PetscInt &ptIdx, const PetscInt &dof, PetscInt &idx) const
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    if ((bIdx < 0) || (bIdx >= nBodies))
        SETERRQ2(comm, PETSC_ERR_ARG_SIZ,
                "Body index %d is out of range. Total number of bodies is %d.",
                bIdx, nBodies);

    ierr = bodies[bIdx]->getGlobalIndex(ptIdx, dof, idx); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // getGlobalIndex


PetscErrorCode BodyPackBase::getGlobalIndex(
        const PetscInt &bIdx, const MatStencil &s, PetscInt &idx) const
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    ierr = getGlobalIndex(bIdx, s.i, s.c, idx); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // getGlobalIndex


PetscErrorCode BodyPackBase::getPackedGlobalIndex(const PetscInt &bIdx, 
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
        idx += bodies[b]->nLclAllProcs[p];

    idx += (unPackedIdx - bodies[bIdx]->offsetsAllProcs[p]);

    PetscFunctionReturn(0);
} // getPackedGlobalIndex


PetscErrorCode BodyPackBase::getPackedGlobalIndex(
        const PetscInt &bIdx, const MatStencil &s, PetscInt &idx) const
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    ierr = getPackedGlobalIndex(bIdx, s.i, s.c, idx); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // getPackedGlobalIndex


PetscErrorCode BodyPackBase::calculateAvgForces(
        const Vec &f, type::RealVec2D &fAvg)
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
        ierr = bodies[i]->calculateAvgForces(unPacked[i], fAvg[i]); CHKERRQ(ierr);
    }

    ierr = DMCompositeRestoreAccessArray(
            dmPack, f, nBodies, nullptr, unPacked.data()); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // calculateAvgForces


PetscErrorCode createBodyPack(const type::Mesh &mesh,
        const YAML::Node &node, type::BodyPack &bodies)
{
    PetscFunctionBeginUser;
    
    bodies = std::make_shared<BodyPackBase>(mesh, node);
    
    PetscFunctionReturn(0);
} // createBodyPack
} // end of namespace body
} // end of namespace petibm
