/**
 * \file singlebodypoints.cpp
 * \brief Implementation of body::SingleBodyPoints.
 * \author Anush Krishnan (anus@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \copyright MIT.
 */


// STL
# include <algorithm>

// PetIBM
# include <petibm/singlebodypoints.h>
# include <petibm/io.h>


namespace petibm
{
namespace body
{


SingleBodyPoints::SingleBodyPoints(const type::Mesh &inMesh, 
        const std::string &inName, const std::string &inFile) :
    SingleBodyBase(inMesh, inName, inFile)
{
    init(inMesh, inName, inFile);
}


PetscErrorCode SingleBodyPoints::init(const type::Mesh &inMesh, 
        const std::string &inName, const std::string &inFile)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    // read coordinates from the given file; nPts and coords are set up here
    ierr = io::readLagrangianPoints(file, nPts, coords); CHKERRQ(ierr);
    
    // check if the dimension of coordinates matches
    if ((unsigned)dim != coords[0].size())
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_READ,
                "The dimension of Lagrangian points are different than that "
                "of the background mesh!\n");

    // create a distributed 1D DMDA with DoF equal to dim; nLclPts, bgPt, edPt,
    // da, nLclAllProcs, and offsetsAllProcs are set up here
    ierr = createDMDA(); CHKERRQ(ierr);

    // initialize meshIdx, which only contains background mesh indices of local
    // Lagrangian points. The indices are defined by pressure cell.
    meshIdx = type::IntVec2D(nLclPts, type::IntVec1D(dim, 0));
    
    // set up background mesh indices for Lagrangian points owned locally
    // meshIdx is defined here.
    ierr = updateMeshIdx(); CHKERRQ(ierr);

    // create info string
    ierr = createInfoString(); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}


PetscErrorCode SingleBodyPoints::createDMDA()
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    DMDALocalInfo       lclInfo;

    ierr = DMDACreate1d(comm, DM_BOUNDARY_NONE, 
            nPts, dim, 0, nullptr, &da); CHKERRQ(ierr);
    
    ierr = DMSetUp(da); CHKERRQ(ierr);

    ierr = DMDAGetLocalInfo(da, &lclInfo); CHKERRQ(ierr);

    // copy necessary local info
    bgPt = lclInfo.xs;
    nLclPts = lclInfo.xm;
    edPt = bgPt + nLclPts;

    // gather local info from other processes
    ierr = MPI_Allgather(&nLclPts, 1, MPIU_INT,
            nLclAllProcs.data(), 1, MPIU_INT, comm); CHKERRQ(ierr);

    // each point has "dim" degree of freedom, so we have to multiply that
    for(auto &it: nLclAllProcs) it *= dim;

    // calculate the offset of the un-packed DM
    for(PetscMPIInt r=mpiSize-1; r>0; r--) offsetsAllProcs[r] = nLclAllProcs[r-1];
    for(PetscMPIInt r=1; r<mpiSize; r++) offsetsAllProcs[r] += offsetsAllProcs[r-1];

    PetscFunctionReturn(0);
}


PetscErrorCode SingleBodyPoints::updateMeshIdx()
{
    PetscFunctionBeginUser;

    // loop through points owned locally and find indices
    for(PetscInt i=bgPt, c=0; i<edPt; ++i, ++c)
    {
        for(PetscInt d=0; d<dim; ++d)
        {
            if (mesh->min[d] >= coords[i][d] || mesh->max[d] <= coords[i][d])
            {
                SETERRQ3(PETSC_COMM_WORLD, PETSC_ERR_MAX_VALUE,
                        "body coordinate %g is outside domain [%g, %g] !",
                        coords[i][d], mesh->min[d], mesh->max[d]);
            }

            meshIdx[c][d] = std::upper_bound(
                    mesh->coord[4][d], mesh->coord[4][d]+mesh->n[4][d],
                    coords[i][d]) - mesh->coord[4][d] - 1;
        }
    }

    PetscFunctionReturn(0);
}


PetscErrorCode SingleBodyPoints::createInfoString()
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

        ss << "\tInput mesh file: " << file << std::endl << std::endl;
        
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

    ierr = MPI_Barrier(comm); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


PetscErrorCode SingleBodyPoints::findProc(const PetscInt &i, PetscMPIInt &p) const
{
    PetscFunctionBeginUser;

    if ((i < 0) || (i >= nPts))
        SETERRQ2(comm, PETSC_ERR_ARG_SIZ, 
                "Index %d of Lagrangian point on the body %s is out of range.",
                i, name.c_str());

    // find the process that own THE 1ST DoF OF THE POINT i
    p = std::upper_bound(offsetsAllProcs.begin(), 
            offsetsAllProcs.end(), i*dim) - offsetsAllProcs.begin() - 1;

    PetscFunctionReturn(0);
}


PetscErrorCode SingleBodyPoints::getGlobalIndex(
        const PetscInt &i, const PetscInt &dof, PetscInt &idx) const
{
    PetscFunctionBeginUser;

    if ((i < 0) || (i >= nPts))
        SETERRQ2(comm, PETSC_ERR_ARG_SIZ, 
                "Index %d of Lagrangian point on the body %s is out of range.",
                i, name.c_str());

    if ((dof < 0) || (dof >= dim))
        SETERRQ2(comm, PETSC_ERR_ARG_SIZ,
                "DoF %d is not correct. The dimension is %d.", dof, dim);

    // for single body DM, the global is simple due to we use 1D DMDA.
    idx = i * dim + dof;
    
    PetscFunctionReturn(0);
}


PetscErrorCode SingleBodyPoints::getGlobalIndex(
        const MatStencil &s, PetscInt &idx) const
{
    PetscFunctionBeginUser;

    PetscErrorCode  ierr;

    ierr = getGlobalIndex(s.i, s.c, idx); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}


PetscErrorCode SingleBodyPoints::calculateAvgForces(
        const Vec &f, type::RealVec1D &fAvg) const
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    PetscReal           **fArry;

    type::RealVec1D    fAvgLocal(dim, 0.0);

    ierr = DMDAVecGetArrayDOF(da, f, &fArry); CHKERRQ(ierr);

    for(PetscInt i=bgPt; i<edPt; ++i)
    {
        for(PetscInt dof=0; dof<dim; ++dof)
        {
            fAvgLocal[dof] -= fArry[i][dof]; // fArray is the force applied to fluid
        }
    }
    ierr = MPI_Barrier(comm); CHKERRQ(ierr);
    
    ierr = MPI_Allreduce(fAvgLocal.data(), fAvg.data(), dim, 
            MPIU_REAL, MPI_SUM, comm); CHKERRQ(ierr);

    ierr = DMDAVecRestoreArrayDOF(da, f, &fArry); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

} // end of namespace body
} // end of namespace petibm
