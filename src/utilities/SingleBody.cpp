/***************************************************************************//**
 * \file SingleBody.cpp
 * \author Anush Krishnan (anus@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Definition of the members of `SingleBody`.
 */


// STL
# include <fstream>
# include <sstream>
# include <algorithm>

// PetIBM
# include "misc.h"
# include "SingleBody.h"


/** \copydoc SingleBody::SingleBody(). */
SingleBody::SingleBody() = default;


/** \copydoc SingleBody::~SingleBody. */
SingleBody::~SingleBody() = default;


/** \copydoc SingleBody::SingleBody(const CartesianMesh &, const std::string &). */
SingleBody::SingleBody(const CartesianMesh &_mesh, 
        const std::string &file, const std::string &_name)
{
    init(_mesh, file, _name);
}


/** \copydoc SingleBody::SingleBogy(const CartesianMesh &, const types::RealVec2D &). */
SingleBody::SingleBody(const CartesianMesh &_mesh, 
        const types::RealVec2D &_coords, const std::string &_name)
{
    init(_mesh, _coords, _name);
}


/** \copydoc SingleBody::init(const CartesianMesh &, const std::string &). */
PetscErrorCode SingleBody::init(const CartesianMesh &_mesh, 
        const std::string &file, const std::string &_name)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    // initialize MPI info and mesh
    ierr = preInit(_mesh, _name); CHKERRQ(ierr);

    // read coordinates from the given file; nPts and coords are set up here
    ierr = readFromFile(file); CHKERRQ(ierr);

    // initialize others
    ierr = postInit(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/** \copydoc SingleBody::init(const CartesianMesh &, const types::RealVec2D &). */
PetscErrorCode SingleBody::init(const CartesianMesh &_mesh, 
        const types::RealVec2D &_coords, const std::string &_name)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    // initialize MPI info and mesh
    ierr = preInit(_mesh, _name); CHKERRQ(ierr);

    // copy cooedinates and set total number of Lagrangian points
    // TODO: check if the input _coords is correct?
    coords = _coords;
    nPts = coords[0].size();

    // initialize others
    ierr = postInit(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/** \copydoc SingleBody::preInit. */
PetscErrorCode SingleBody::preInit(
        const CartesianMesh &_mesh, const std::string &_name)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    // set up the name
    name = _name;

    // save the address of input mesh instance
    // note: a bad practice for shared_ptr
    mesh = std::shared_ptr<const CartesianMesh>(&_mesh, [](const CartesianMesh*){});

    // save MPI information from the mesh
    comm = mesh->comm;
    mpiSize = mesh->mpiSize;
    mpiRank = mesh->mpiRank;

    // set the dimension
    dim = mesh->dim;

    // allocate vectors
    nLclAllProcs = types::IntVec1D(mpiSize, 0);
    offsetsAllProcs = types::IntVec1D(mpiSize, 0);

    PetscFunctionReturn(0);
}


/** \copydoc SingleBody::postInit. */
PetscErrorCode SingleBody::postInit()
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    // create a distributed 1D DMDA with DoF equal to dim; nLclPts, bgPt, edPt,
    // da, nLclAllProcs, and offsetsAllProcs are set up here
    ierr = createDMDA(); CHKERRQ(ierr);

    // set up background mesh indices for Lagrangian points owned locally
    // meshIdx is defined here.
    ierr = findCellIdx(); CHKERRQ(ierr);

    // create info string
    ierr = createInfoString(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/** \copydoc SingleBody::readFromFile. */
PetscErrorCode SingleBody::readFromFile(const std::string &file)
{
    PetscFunctionBeginUser;

    std::ifstream   inFile(file);
    std::string     line;
    PetscInt        c = 0;

    // check if we successfully open the file
    if (! inFile.good()) SETERRQ1(*comm, PETSC_ERR_FILE_READ,
            "Opening or reading body file %s failed!", file.c_str());


    // read the total number of points
    if (std::getline(inFile, line)) 
    {
        std::stringstream   sline(line);

        if (! (sline >> nPts))
            SETERRQ1(*comm, PETSC_ERR_FILE_READ, 
                    "Can't read the total number of points in file %s !\n",
                    file.c_str());

        if (sline.peek() != EOF) 
            SETERRQ1(*comm, PETSC_ERR_FILE_READ, 
                    "The first line in file %s contains more than one integer. "
                    "Please check the format.\n", file.c_str());
    }
    else
        SETERRQ1(*comm, PETSC_ERR_FILE_READ,
                "Error while reading the first line in file %s !\n", file.c_str());


    // initialize the size of coordinate array
    coords = types::RealVec2D(3, types::RealVec1D(nPts, 0.0));

    // loop through all points and save their coordinates
    while(std::getline(inFile, line))
    {
        std::stringstream   sline(line);

        for(PetscInt d=0; d<dim; ++d)
            if (! (sline >> coords[d][c]))
                SETERRQ2(*comm, PETSC_ERR_FILE_READ,
                        "The number of doubles at line %d in file %s does not "
                        "match the dimension.\n", c+2, file.c_str());

        if (sline.peek() != EOF)
            SETERRQ2(*comm, PETSC_ERR_FILE_READ,
                        "The number of doubles at line %d in file %s does not "
                        "match the dimension.\n", c+2, file.c_str());

        c += 1;
    }

    // close the file
    inFile.close();

    // check the total number of points read
    if (c != nPts)
        SETERRQ1(*comm, PETSC_ERR_FILE_READ,
                "The total number of coordinates read in does not match the "
                "number specified at the first line in file %s !\n",
                file.c_str());

    PetscFunctionReturn(0);
}


/** \copydoc SingleBody::findCellIdx. */
PetscErrorCode SingleBody::findCellIdx()
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    // initialize meshIdx, which only contains background mesh indices of local
    // Lagrangian points. The indices are defined by pressure cell.
    meshIdx = types::IntVec2D(3, types::IntVec1D(nLclPts, 0));

    // loop through points owned locally and find indices
    for(PetscInt d=0; d<dim; ++d)
    {
        PetscInt    c = 0;
        for(PetscInt i=bgPt; i<edPt; ++i)
        {
            ierr = misc::findCell1D(coords[d][i], 
                    mesh->coord[4][d], meshIdx[d][c]); CHKERRQ(ierr);

            c += 1;
        }
    }

    PetscFunctionReturn(0);
}


/** \copydoc SingleBody::createDMDA. */
PetscErrorCode SingleBody::createDMDA()
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    DMDALocalInfo       lclInfo;

    ierr = DMDACreate1d(*comm, DM_BOUNDARY_NONE, 
            nPts, dim, 0, nullptr, &da); CHKERRQ(ierr);

    ierr = DMDAGetLocalInfo(da, &lclInfo); CHKERRQ(ierr);

    // copy necessary local info
    bgPt = lclInfo.xs;
    nLclPts = lclInfo.xm;
    edPt = bgPt + nLclPts;

    // gather local info from other processes
    ierr = MPI_Allgather(&nLclPts, 1, MPIU_INT,
            nLclAllProcs.data(), 1, MPIU_INT, *comm); CHKERRQ(ierr);

    // each point has dim degree of freedom, so we have to multiply that
    for(auto &it: nLclAllProcs) it *= dim;

    // calculate the offset of the un-packed DM
    for(PetscMPIInt r=mpiSize-1; r>0; r--) offsetsAllProcs[r] = nLclAllProcs[r-1];
    for(PetscMPIInt r=1; r<mpiSize; r++) offsetsAllProcs[r] += offsetsAllProcs[r-1];

    PetscFunctionReturn(0);
}


/** \copydoc SingleBody::createInfoString. */
PetscErrorCode SingleBody::createInfoString()
{
    PetscFunctionBeginUser;

    PetscErrorCode          ierr;

    std::stringstream       ss;

    // only rank 0 prepares the header of info string
    if (mpiRank == 0)
    {

        ss << std::string(80, '=') << std::endl;
        ss << "Body " << name << ":" << std::endl;
        ss << std::string(80, '=') << std::endl;

        ss << "\tDimension: " << dim << std::endl << std::endl;

        ss << "\tTotal number of Lagrangian points: "
            << nPts << std::endl << std::endl;

        ss << "\tBody is distributed to " << mpiSize
            << " processes" << std::endl << std::endl;

        ss << "\tDistribution of Lagrangian points:" << std::endl << std::endl;
    }

    ss << "\t\tRank " << mpiRank << ":" << std::endl;
    ss << "\t\t\tNumber of points: " << nLclPts << std::endl;
    ss << "\t\t\tRange of points: [" << bgPt << ", " << edPt << ")" << std::endl;

    info = ss.str();

    ierr = MPI_Barrier(*comm); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/** \copydoc SingleBody::printInfo. */
PetscErrorCode SingleBody::printInfo()
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    ierr = PetscSynchronizedPrintf(*comm, info.c_str()); CHKERRQ(ierr);
    ierr = PetscSynchronizedFlush(*comm, PETSC_STDOUT); CHKERRQ(ierr);
    ierr = PetscPrintf(*comm, "\n"); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/** \copydoc SingleBody::findProc. */
PetscErrorCode SingleBody::findProc(const PetscInt &i, PetscMPIInt &p)
{
    PetscFunctionBeginUser;

    if ((i < 0) || (i >= nPts))
        SETERRQ2(*comm, PETSC_ERR_ARG_SIZ, 
                "Index %d of Lagrangian point on the body %s is out of range.",
                i, name.c_str());

    // find the process that own THE 1ST DoF OF THE POINT i
    p = std::upper_bound(offsetsAllProcs.begin(), 
            offsetsAllProcs.end(), (i-1)*dim+1) - offsetsAllProcs.begin() - 1;

    PetscFunctionReturn(0);
}


/** \copydoc SingleBody::findGlobalIndex(
 *           const PetscInt &, const PetscInt &, PetscInt &). */
PetscErrorCode SingleBody::getGlobalIndex(
        const PetscInt &i, const PetscInt &dof, PetscInt &idx)
{
    PetscFunctionBeginUser;

    if ((i < 0) || (i >= nPts))
        SETERRQ2(*comm, PETSC_ERR_ARG_SIZ, 
                "Index %d of Lagrangian point on the body %s is out of range.",
                i, name.c_str());

    if ((dof < 0) || (dof >= dim))
        SETERRQ2(*comm, PETSC_ERR_ARG_SIZ,
                "DoF %d is not correct. The dimension is %d.", dof, dim);

    // for single body DM, the global is simple due to we use 1D DMDA.
    idx = (i - 1) * dim + dof;
    
    PetscFunctionReturn(0);
}


/** \copydoc SingleBody::findGlobalIndex(const MatStencil &s, PetscInt &idx). */
PetscErrorCode SingleBody::getGlobalIndex(const MatStencil &s, PetscInt &idx)
{
    PetscFunctionBeginUser;

    PetscErrorCode  ierr;

    ierr = getGlobalIndex(s.i, s.c, idx); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}
