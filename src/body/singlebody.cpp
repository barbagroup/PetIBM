/**
 * \file singlebody.cpp
 * \brief Implementations of body::SingleBodyBase and factory function.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 */

// PetIBM
#include <petibm/io.h>
#include <petibm/singlebody.h>
#include <petibm/singlebodypoints.h>

namespace petibm
{
namespace body
{
SingleBodyBase::SingleBodyBase(const type::Mesh &inMesh,
                               const std::string &inName,
                               const std::string &inFile)
{
    // set up the name
    name = inName;

    // store the path of input file
    file = inFile;

    // save reference to the background mesh
    mesh = inMesh;

    // save MPI information from the mesh
    comm = mesh->comm;
    mpiSize = mesh->mpiSize;
    mpiRank = mesh->mpiRank;

    // save the dimension
    dim = mesh->dim;

    // allocate vectors
    nLclAllProcs = type::IntVec1D(mpiSize, 0);
    offsetsAllProcs = type::IntVec1D(mpiSize, 0);
}  // SingleBodyBase

SingleBodyBase::~SingleBodyBase()
{
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    PetscBool finalized;

    ierr = PetscFinalized(&finalized); CHKERRV(ierr);
    if (finalized) return;

    ierr = DMDestroy(&da); CHKERRV(ierr);
    comm = MPI_COMM_NULL;
}  // ~SingleBodyBase

PetscErrorCode SingleBodyBase::destroy()
{
    PetscFunctionBeginUser;
    PetscErrorCode ierr;

    dim = -1;
    name = file = info = "";
    nPts = nLclPts = bgPt = edPt = 0;
    type::RealVec2D().swap(coords);
    type::IntVec2D().swap(meshIdx);
    ierr = DMDestroy(&da); CHKERRQ(ierr);
    comm = MPI_COMM_NULL;
    mpiSize = mpiRank = 0;
    mesh.reset();
    type::IntVec1D().swap(nLclAllProcs);
    type::IntVec1D().swap(offsetsAllProcs);

    PetscFunctionReturn(0);
}  // destroy

PetscErrorCode SingleBodyBase::printInfo() const
{
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    ierr = io::print(info); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}  // printInfo

PetscErrorCode createSingleBody(const type::Mesh &mesh, const std::string &type,
                                const std::string &name,
                                const std::string &file, type::SingleBody &body)
{
    PetscFunctionBeginUser;

    if (type == "points")
        body = std::make_shared<SingleBodyPoints>(mesh, name, file);
    else
        SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                 "The type of mesh file \"%s\" is not recognized!\n",
                 type.c_str());

    PetscFunctionReturn(0);
}  // createSingleBody

}  // end of namespace body
}  // end of namespace petibm
