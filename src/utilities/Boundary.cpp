/***************************************************************************//**
 * \file Boundary.h
 * \author Anush Krishnan (anus@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Definition of the class `Boundary`.
 */


// here goes headers from our PetIBM
#include "Boundary.h"


namespace petibm
{
namespace utilities
{

using namespace types;

/** \copydoc Boundary::Boundary() */
Boundary::Boundary() = default;


/** \copydoc Boundary::~Boundary() */
Boundary::~Boundary() = default;


/** \copydoc Boundary::Boundary(const CartesianMesh &) */
Boundary::Boundary(const CartesianMesh &mesh) { init(mesh); }


/** \copydoc Boundary::init(const Cartesian &) */
PetscErrorCode Boundary::init(const CartesianMesh &_mesh)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    // create a shared pointer to mesh; bad practice...
    mesh = std::shared_ptr<const CartesianMesh>(&_mesh, [](const CartesianMesh*){}); 

    // obtain MPI information from CartesianMesh object
    comm = mesh->comm;
    mpiSize = mesh->mpiSize;
    mpiRank = mesh->mpiRank;

    // set dim
    dim = mesh->dim;

    // initialize the size of bds (Question: is barrier necessary?)
    bds.resize(dim * 2);

    bds[0].init(*mesh, BCLoc::XMINUS);
    bds[1].init(*mesh, BCLoc::XPLUS);
    bds[2].init(*mesh, BCLoc::YMINUS);
    bds[3].init(*mesh, BCLoc::YPLUS);

    if (dim == 3)
    {
        bds[4].init(*mesh, BCLoc::ZMINUS);
        bds[5].init(*mesh, BCLoc::ZPLUS);
    }

    PetscFunctionReturn(0);
}


/** \copydoc Boundary::setGhostICs(const Solutions &). */
PetscErrorCode Boundary::setGhostICs(const Solutions &soln)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    for(auto &bd: bds) { ierr = bd.setGhostICs(soln); CHKERRQ(ierr); }

    ierr = MPI_Barrier(*comm); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/** \copydoc Boundary::updateCoeffs */
PetscErrorCode Boundary::updateEqs(const Solutions &soln, const PetscReal &dt)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    for(auto &it: bds) { ierr = it.updateEqs(soln, dt); CHKERRQ(ierr); }

    ierr = MPI_Barrier(*comm); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/** \copydoc Boundary::updateGhosts */
PetscErrorCode Boundary::updateGhostValues(const Solutions &soln)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    for(auto &it: bds) { ierr = it.updateGhostValues(soln); CHKERRQ(ierr); }

    ierr = MPI_Barrier(*comm); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/** \copydoc Boundary::copyValues2LocalVecs */
PetscErrorCode Boundary::copyValues2LocalVecs(std::vector<Vec> &lclVecs) const
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    for(auto &it: bds) { ierr = it.copyValues2LocalVecs(lclVecs); CHKERRQ(ierr); }

    ierr = MPI_Barrier(*comm); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

} // end of namespace utilities
} // end of namespace petibm
