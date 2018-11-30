/**
 * \file singlebody.cpp
 * \brief Implementations of body::SingleBodyBase and factory function.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 */

#include <petibm/io.h>
#include <petibm/singlebody.h>
#include <petibm/singlebodypoints.h>

namespace petibm
{
namespace body
{
SingleBodyBase::SingleBodyBase(const MPI_Comm &inComm, const PetscInt &inDim,
                               const std::string &inName,
                               const std::string &inFilePath)
{
    // set up the name
    name = inName;

    // store the path of input file
    filePath = inFilePath;

    // store MPI information
    comm = inComm;
    MPI_Comm_size(comm, &mpiSize);
    MPI_Comm_rank(comm, &mpiRank);

    // save the dimension
    dim = inDim;

    // allocate vectors
    nLclAllProcs = type::IntVec1D(mpiSize, 0);
    offsetsAllProcs = type::IntVec1D(mpiSize, 0);
}  // SingleBodyBase

SingleBodyBase::~SingleBodyBase()
{
    PetscErrorCode ierr;
    PetscBool finalized;

    PetscFunctionBeginUser;

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
    name = filePath = info = "";
    nPts = nLclPts = bgPt = edPt = 0;
    type::RealVec2D().swap(coords);
    type::IntVec2D().swap(meshIdx);
    ierr = DMDestroy(&da); CHKERRQ(ierr);
    comm = MPI_COMM_NULL;
    mpiSize = mpiRank = 0;
    type::IntVec1D().swap(nLclAllProcs);
    type::IntVec1D().swap(offsetsAllProcs);

    PetscFunctionReturn(0);
}  // destroy

PetscErrorCode SingleBodyBase::printInfo() const
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = io::print(info); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // printInfo

PetscErrorCode createSingleBody(const MPI_Comm &comm, const PetscInt &dim,
                                const std::string &type,
                                const std::string &name,
                                const std::string &filePath,
                                type::SingleBody &body)
{
    PetscFunctionBeginUser;

    if (type == "points")
        body = std::make_shared<SingleBodyPoints>(comm, dim, name, filePath);
    else
        SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                 "The type of body file \"%s\" is not recognized!\n",
                 type.c_str());

    PetscFunctionReturn(0);
}  // createSingleBody

}  // end of namespace body

}  // end of namespace petibm
