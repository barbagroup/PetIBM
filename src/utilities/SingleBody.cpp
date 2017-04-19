/***************************************************************************//**
 * \file SingleBody.cpp
 * \author Anush Krishnan (anus@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Definition of the members of `SingleBody`.
 */


// STL
# include <fstream>

// PetIBM
# include "misc.h"
# include "SingleBody.h"


/** \copydoc SingleBody::SingleBody(). */
SingleBody::SingleBody() = default;


/** \copydoc SingleBody::~SingleBody. */
SingleBody::~SingleBody() = default;


/** \copydoc SingleBody::SingleBody(const CartesianMesh &, const std::string &). */
SingleBody::SingleBody(
        const CartesianMesh &_mesh, const std::string &file)
{
    init(_mesh, file);
}


/** \copydoc SingleBody::SingleBogy(const CartesianMesh &, const types::RealVec2D &). */
SingleBody::SingleBody(
        const CartesianMesh &_mesh, const types::RealVec2D &_coords)
{
    init(_mesh, _coords);
}


/** \copydoc SingleBody::init(const CartesianMesh &, const std::string &). */
PetscErrorCode SingleBody::init(
        const CartesianMesh &_mesh, const std::string &file)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    // initialize MPI info and mesh
    ierr = preInit(_mesh); CHKERRQ(ierr);

    // read coordinates from the given file; nPts and coords are set up here
    ierr = readFromFile(file); CHKERRQ(ierr);

    // initialize others
    ierr = postInit(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/** \copydoc SingleBody::init(const CartesianMesh &, const types::RealVec2D &). */
PetscErrorCode SingleBody::init(
        const CartesianMesh &_mesh, const types::RealVec2D &_coords)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    // initialize MPI info and mesh
    ierr = preInit(_mesh); CHKERRQ(ierr);

    // copy cooedinates and set total number of Lagrangian points
    // TODO: check if the input _coords is correct?
    coords = _coords;
    nPts = coords[0].size();

    // initialize others
    ierr = postInit(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/** \copydoc SingleBody::preInit. */
PetscErrorCode SingleBody::preInit(const CartesianMesh &_mesh)
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

    PetscFunctionReturn(0);
}


/** \copydoc SingleBody::postInit. */
PetscErrorCode SingleBody::postInit()
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    // create a distributed 1D DMDA with DoF equal to dim; nLclPts, bg, ed and
    // da are set up here
    ierr = createDMDA(); CHKERRQ(ierr);

    // set up background mesh indices for Lagrangian points owned locally
    // meshIdx is defined here.
    ierr = findCellIdx(); CHKERRQ(ierr);

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
        for(PetscInt i=bg; i<ed; ++i)
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
    bg = lclInfo.xs;
    nLclPts = lclInfo.xm;
    ed = bg + nLclPts;

    PetscFunctionReturn(0);
}
