/***************************************************************************//**
 * \file singleboundary.cpp
 * \author Anush Krishnan (anus@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Definition of the member functions of class `SingleBoundary`.
 */

// here goes headers from our PetIBM
# include <petibm/misc.h>
# include <petibm/singleboundary.h>


namespace petibm
{
namespace boundary
{

SingleBoundaryBase::SingleBoundaryBase(const type::Mesh &inMesh,
        const type::BCLoc &bcLoc, const type::Field &inField, const PetscReal &inValue)
{
    init(inMesh, bcLoc, inField, inValue);
}


PetscErrorCode SingleBoundaryBase::init(const type::Mesh &inMesh,
        const type::BCLoc &bcLoc, const type::Field &inField, const PetscReal &inValue)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    // create a shared pointer to mesh; bad practice...
    mesh = inMesh; 

    // obtain MPI information from CartesianMesh object
    comm = mesh->comm;
    mpiSize = mesh->mpiSize;
    mpiRank = mesh->mpiRank;

    // set dim
    dim = mesh->dim;

    // set the location
    loc = bcLoc;
    
    // set field
    field = inField;
    
    // set value
    value = inValue;
    
    // set normal
    normal = ((int(loc)%2) == 0) ? -1.0 : 1.0;

    // set onThisProc
    ierr = misc::checkBoundaryProc(mesh->da[int(field)],
            mesh->n[int(field)], loc, onThisProc); CHKERRQ(ierr);

    // for processes on this boundary, set up IDs and values
    if (onThisProc)
    {
        ierr = misc::getGhostPointList(mesh, field, loc, points); CHKERRQ(ierr);
    }

    ierr = MPI_Barrier(comm); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


PetscErrorCode SingleBoundaryBase::setGhostICs(const Vec &vec)
{
    PetscFunctionBeginUser;

    PetscErrorCode  ierr;
    
    PetscReal       targetValue;
    
    for(auto &it: points)
    {
        ierr = VecGetValues(vec, 1,
                &(it.second.targetPackedId), &targetValue); CHKERRQ(ierr);
        
        ierr = setGhostICsKernel(targetValue, it.second); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}


PetscErrorCode SingleBoundaryBase::updateEqs(const Vec &vec, const PetscReal &dt)
{
    PetscFunctionBeginUser;

    PetscErrorCode  ierr;
    
    PetscReal       targetValue;
    
    for(auto &it: points)
    {
        ierr = VecGetValues(vec, 1,
                &(it.second.targetPackedId), &targetValue); CHKERRQ(ierr);
        
        ierr = updateEqsKernel(targetValue, dt, it.second); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}


PetscErrorCode SingleBoundaryBase::updateGhostValues(const Vec &vec)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    PetscReal           targetValue;

    for(auto &it: points)
    {
        ierr = VecGetValues(vec, 1, 
                &(it.second.targetPackedId), &targetValue); CHKERRQ(ierr);

        it.second.value = it.second.a0 * targetValue + it.second.a1;
    }

    PetscFunctionReturn(0);
}


PetscErrorCode SingleBoundaryBase::copyValues2LocalVec(Vec &lclVec)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    for(auto &it: points)
    {
        ierr = VecSetValue(lclVec, it.second.lclId, 
                it.second.value, INSERT_VALUES); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}

} // end of namespace boundary
} // end of namespace petibm
