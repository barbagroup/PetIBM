/**
 * \file boundarysimple.cpp
 * \brief Implementation of boundary::BoundarySimple
 * \author Anush Krishnan (anus@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \copyright MIT.
 */


// here goes headers from our PetIBM
# include <petibm/boundarysimple.h>
# include <petibm/parser.h>


namespace petibm
{
namespace boundary
{
using namespace type;


// constructor
BoundarySimple::BoundarySimple(const type::Mesh &inMesh, const YAML::Node &node)
{
    init(inMesh, node);
}


// underlying initialization function
PetscErrorCode BoundarySimple::init(
        const type::Mesh &inMesh, const YAML::Node &node)
{
    PetscFunctionBeginUser;
    
    PetscErrorCode  ierr;

    // make a shared pointer pointing to associated mesh
    mesh = inMesh; 

    // obtain MPI information from CartesianMesh object
    comm = mesh->comm;
    mpiSize = mesh->mpiSize;
    mpiRank = mesh->mpiRank;

    // set dim
    dim = mesh->dim;

    // initialize the size of bds (Question: is barrier necessary?)
    bds.resize(dim);
    
    type::IntVec2D      bcTypes;
    type::RealVec2D     bcValues;
    ierr = parser::parseBCs(node, bcTypes, bcValues); CHKERRQ(ierr);

    for(PetscInt f=0; f<dim; ++f)
    {
        bds[f].resize(dim*2);
        
        for(PetscInt b=0; b<dim*2; ++b)
        {
            ierr = createSingleBoundary(
                    mesh, type::BCLoc(b), type::Field(f), 
                    bcValues[f][b], type::BCType(bcTypes[f][b]),
                    bds[f][b]); CHKERRQ(ierr);
        }
    }

    PetscFunctionReturn(0);
}


// set intial values to ghost points
PetscErrorCode BoundarySimple::setGhostICs(const type::Solution &soln)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    for(auto &fbd: bds)
    { 
        for(auto &bd: fbd)
        {
            ierr = bd->setGhostICs(soln->UGlobal); CHKERRQ(ierr);
        }
    }

    ierr = MPI_Barrier(comm); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


// update the equations between boundary and ghost points
PetscErrorCode BoundarySimple::updateEqs(
        const type::Solution &soln, const PetscReal &dt)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    for(auto &fbd: bds)
    { 
        for(auto &bd: fbd)
        {
            ierr = bd->updateEqs(soln->UGlobal, dt); CHKERRQ(ierr);
        }
    }

    ierr = MPI_Barrier(comm); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


// update values of ghost points
PetscErrorCode BoundarySimple::updateGhostValues(const type::Solution &soln)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    for(auto &fbd: bds)
    { 
        for(auto &bd: fbd)
        {
            ierr = bd->updateGhostValues(soln->UGlobal); CHKERRQ(ierr);
        }
    }

    ierr = MPI_Barrier(comm); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


// copy values of ghost points to local PETSc Vecs
PetscErrorCode BoundarySimple::copyValues2LocalVecs(std::vector<Vec> &lclVecs) const
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    for(PetscInt f=0; f<dim; ++f)
    { 
        for(auto &bd: bds[f])
        {
            ierr = bd->copyValues2LocalVec(lclVecs[f]); CHKERRQ(ierr);
        }
    }

    ierr = MPI_Barrier(comm); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

} // end of namespace boundary
} // end of namespace petibm
