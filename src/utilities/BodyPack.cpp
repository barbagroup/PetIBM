/***************************************************************************//**
 * \file BodyPack.cpp
 * \author Anush Krishnan (anus@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Definition of the members of class `BodyPack`.
 */


// PetIBM
# include "BodyPack.h"


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
    }

    // create dmPack
    ierr = createDmPack(); CHKERRQ(ierr);

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


/** \copydoc BodyPack::createInfoString. */
PetscErrorCode BodyPack::createInfoString()
{
    PetscFunctionBeginUser;

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

    PetscFunctionReturn(0);
}


/** \copydoc BodyPack::printInfo. */
PetscErrorCode BodyPack::printInfo()
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    ierr = PetscPrintf(*comm, info.c_str()); CHKERRQ(ierr);

    for(PetscInt i=0; i<nBodies; ++i)
    {
        ierr = bodies[i].printInfo(); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}
