/***************************************************************************//**
 * \file cartesianmesh.cpp
 * \author Anush Krishnan (anush@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Implementation of the methods of the class `CartesianMesh`.
 */


// here goes C++ STL
# include <iostream>
# include <sstream>
# include <algorithm>
# include <numeric>

// here goes PETSc headers
# include <petscvec.h>

// here goes headers from our PetIBM
# include <petibm/cartesianmesh.h>
# include <petibm/misc.h>
# include <petibm/parser.h>


namespace petibm
{
namespace mesh
{

using namespace type;


// constructor
CartesianMesh::CartesianMesh(const MPI_Comm &world, const YAML::Node &node)
{
    init(world, node);
}


// default destructor
CartesianMesh::~CartesianMesh() = default;


// initialization
PetscErrorCode CartesianMesh::init(const MPI_Comm &world, const YAML::Node &node)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    // store the address of the communicator
    // note: this is a bad practice; shared_ptr is not for stack variables!!
    comm = world;

    // set rank and size
    ierr = MPI_Comm_size(comm, &mpiSize); CHKERRQ(ierr);
    ierr = MPI_Comm_rank(comm, &mpiRank); CHKERRQ(ierr);


    // set the default sizes and values for members
    // The default values are chosen for 2D cases. We want to eliminate the use
    // of template and also `if` conditions that determine the dimension. So
    // with carefully choosed default values, we can take care of 2D case with 
    // 3D code in many places.
    min = RealVec1D(3, 0.0);
    max = RealVec1D(3, 1.0);
    n = IntVec2D(5, IntVec1D(3, 1));
    coordTrue = RealVec3D(5, RealVec2D(3, RealVec1D(1, 0.0)));
    coord = GhostedVec3D(5, GhostedVec2D(3, nullptr));
    dLTrue = RealVec3D(5, RealVec2D(3, RealVec1D(1, 1.0)));
    dL = GhostedVec3D(5, GhostedVec2D(3, nullptr));
    da = std::vector<DM>(5, PETSC_NULL);
    nProc = IntVec1D(3, PETSC_DECIDE);
    bg = IntVec2D(5, IntVec1D(3, 0));
    ed = IntVec2D(5, IntVec1D(3, 1));
    m = IntVec2D(5, IntVec1D(3, 0));
    ao = std::vector<AO>(4);
    UNLocalAllProcs = IntVec2D(3, IntVec1D(mpiSize, 0));
    UPackNLocalAllProcs = IntVec1D(mpiSize, 0);
    offsetsAllProcs = IntVec2D(3, IntVec1D(mpiSize, 0));
    offsetsPackAllProcs = IntVec1D(mpiSize, 0);


    // index 3 represent pressure mesh; min & max always represent pressure mesh
    ierr = parser::parseMesh(node["mesh"], dim, min, max, n[3], dLTrue[3]); CHKERRQ(ierr);

    // check periodic BC
    IntVec2D          bcTypes;
    RealVec2D         bcValues;
    BoolVec2D         periodic;
    ierr = parser::parseBCs(node, bcTypes, bcValues); CHKERRQ(ierr);
    ierr = misc::checkPeriodicBC(bcTypes, periodic); CHKERRQ(ierr);
    
    // create raw grid information
    ierr = createPressureMesh(); CHKERRQ(ierr);
    ierr = createVertexMesh(); CHKERRQ(ierr);
    ierr = createVelocityMesh(periodic); CHKERRQ(ierr);
    ierr = MPI_Barrier(comm); CHKERRQ(ierr);

    // create PETSc DMs
    ierr = initDMDA(periodic); CHKERRQ(ierr);

    // setup the format of the file that the `write` function will use
    //ierr = setOutputFormat(type); CHKERRQ(ierr);

    // create a std::string that can be used in `printInfo` and output stream
    ierr = createInfoString(); CHKERRQ(ierr);

    // all processes should be syncrinized
    ierr = MPI_Barrier(comm); CHKERRQ(ierr);


    PetscFunctionReturn(0);
}


// create pressure mesh
PetscErrorCode CartesianMesh::createPressureMesh()
{
    PetscFunctionBeginUser;

    // loop through all axes to get the coordinates of pressure points
    for(int i=0; i<dim; ++i)
    {
        coordTrue[3][i].resize(n[3][i]);

        // the following three lines calculate the coord of pressure points
        std::partial_sum(dLTrue[3][i].begin(), 
                dLTrue[3][i].end(), coordTrue[3][i].begin());

        auto f = [i, this] (double &a, double &b) -> double 
                    { return a + this->min[i] - 0.5 * b; };

        std::transform(coordTrue[3][i].begin(), coordTrue[3][i].end(), 
                dLTrue[3][i].begin(), coordTrue[3][i].begin(), f);

        // no ghost points in pressure grid, so dL[0] = dLTrue[0]
        dL[3][i] = &dLTrue[3][i][0];

        // same for coordinates
        coord[3][i] = &coordTrue[3][i][0];
    }

    // total number of pressure points
    pN = n[3][0] * n[3][1] * n[3][2];   // for 2D, n[3][2] should be 1

    // in 2D simulation, we still use the varaible dL in z-direction
    if (dim == 2)
    {
        dL[3][2] = &dLTrue[3][2][0];
        coord[3][2] = &coordTrue[3][2][0];
    }

    PetscFunctionReturn(0);
}


// create vertices information
PetscErrorCode CartesianMesh::createVertexMesh()
{
    PetscFunctionBeginUser;

    // loop through all axes to get the coordinates of mesh vertexes
    // note: index 3 means pressure mesh; 4 means vertexes
    for(int i=0; i<dim; ++i)
    {
        // number of vertexes is one more than pressure cells
        n[4][i] = n[3][i] + 1;
        coordTrue[4][i].resize(n[4][i]);

        // create coordinates of vertexes
        std::partial_sum(dLTrue[3][i].begin(), dLTrue[3][i].end(), 
                coordTrue[4][i].begin()+1); // first assume coord[4][i][0]=0.0

        // shift aacording to real coord[4][i][0], which is min
        std::for_each(coordTrue[4][i].begin(), coordTrue[4][i].end(),
                [i, this] (double &a) -> void { a += this->min[i]; });

        // no ghost points in vertex mesh, so coord[4][i][0]=coordTrue[4][i][0]
        coord[4][i] = &coordTrue[4][i][0];
    }

    // in 2D simulation, we still use the varaible coord in z-direction
    if (dim == 2) coord[4][2] = &coordTrue[4][2][0];

    PetscFunctionReturn(0);
}


// create velocity mesh
PetscErrorCode CartesianMesh::createVelocityMesh(const BoolVec2D &periodic)
{
    PetscFunctionBeginUser;

    // initialize UN
    UN = 0;

    // loop through all velocity component
    // note: index 0, 1, 2 represent u, v, w respectively
    for(int comp=0; comp<dim; comp++)
    {
        // loop through each direction, i.e., x, y, z directions
        for(int dir=0; dir<dim; ++dir)
        {
            // when the direction corresponding to velocity component
            if (dir == comp)
            {
                // there will no velocity point on boundary (it's a ghost point)
                n[comp][dir] = n[3][dir] - 1;

                // we store the dL and coord of ghost points, so the size is n+2
                dLTrue[comp][dir].resize(n[comp][dir]+2);
                coordTrue[comp][dir].resize(n[comp][dir]+2);

                // coordinates will match the those of the vertex grid
                coordTrue[comp][dir] = coordTrue[4][dir];

                // cell size will be a half of the sum of adjacent pressure cell
                auto f = [] (const PetscReal &x, const PetscReal &y)
                    -> PetscReal { return 0.5 * (x + y); };

                // note: adjacent can not automatically ignore the first element
                // note: pressure doesn't have ghost points, while velocity does.
                std::adjacent_difference(
                        dLTrue[3][dir].begin(), dLTrue[3][dir].end(),
                        dLTrue[comp][dir].begin(), f);


                // if BC in this direction is periodic, we append one more 
                // velocity point at the end
                if (periodic[comp][dir])
                {
                    // add 1 to the number of valid grid point
                    n[comp][dir] += 1;

                    // the space for right ghost point in dLTrue is now used to 
                    // store the grid point on periodic BC
                    dLTrue[comp][dir].back() = 
                        f(dLTrue[3][dir][0], dLTrue[3][dir].back());

                    // the smaller ghost point will be the same as the last
                    // valid point
                    dLTrue[comp][dir][0] = dLTrue[comp][dir].back();

                    // we have to create one more space for the larger ghost
                    // point, and its dL is the same as the 1st valid point
                    dLTrue[comp][dir].push_back(dLTrue[comp][dir][1]);

                    // add the coordinate of the extra point, which is the first
                    // interior point add a shift of domain size.
                    coordTrue[comp][dir].push_back(
                            max[dir] + dLTrue[3][dir].front());
                }
                // for other cases, we simply get the dL for the larger ghost
                // note: the smaller ghost point has already got value from
                // adjacent_difference, which is the dL of the 1st pressure cell
                else
                {
                    dLTrue[comp][dir].back() = dLTrue[3][dir].back();
                }
            }
            // other directions
            else
            {
                n[comp][dir] = n[3][dir];

                // we store the dL of ghost points, so the size is n+2
                dLTrue[comp][dir].resize(n[comp][dir]+2);
                coordTrue[comp][dir].resize(n[comp][dir]+2);

                // coordinate and cel size will match that of pressure point,
                // except the two ghost points
                std::copy(coordTrue[3][dir].begin(), coordTrue[3][dir].end(),
                        coordTrue[comp][dir].begin()+1);
                std::copy(dLTrue[3][dir].begin(), dLTrue[3][dir].end(), 
                        dLTrue[comp][dir].begin()+1);

                // get the dL for the two ghost points
                // for periodic BCs, it's the same as the 1st cell on the other side
                if (periodic[comp][dir])
                {
                    // set the coordinate of the two ghost points
                    coordTrue[comp][dir].front() = 
                        min[dir] - dLTrue[3][dir].back() / 2.0;
                    coordTrue[comp][dir].back() = 
                        max[dir] + dLTrue[3][dir].front() / 2.0;

                    // dL
                    dLTrue[comp][dir].front() = dLTrue[3][dir].back();
                    dLTrue[comp][dir].back() = dLTrue[3][dir].front();
                }
                else // for other BCs, it's just a mirror of the first/last cells
                {
                    // set the coordinate of the two ghost points
                    coordTrue[comp][dir].front() = 
                        min[dir] - dLTrue[3][dir].front() / 2.0;
                    coordTrue[comp][dir].back() = 
                        max[dir] + dLTrue[3][dir].back() / 2.0;

                    // dL
                    dLTrue[comp][dir].front() = dLTrue[3][dir].front();
                    dLTrue[comp][dir].back() = dLTrue[3][dir].back();
                }
            }


            // dL and coord will point to the 1st valid grid point, so we can 
            // access the 1st ghost point with the index -1, e.g., dL[0][0][-1]
            dL[comp][dir] = &dLTrue[comp][dir][1];
            coord[comp][dir] = &coordTrue[comp][dir][1];
        }

        // add to UN; for 2D, n in z-direction should be 1
        UN += (n[comp][0] * n[comp][1] * n[comp][2]); 
    }


    // for 2D simulation, we still use dL in z-direction in our code
    if (dim == 2)
    {
        dL[0][2] = &dLTrue[0][2][0]; dL[1][2] = &dLTrue[1][2][0];
        dL[2][0] = &dLTrue[2][0][0]; dL[2][1] = &dLTrue[2][1][0];
        dL[2][2] = &dLTrue[2][2][0];

        coord[0][2] = &coordTrue[0][2][0]; coord[1][2] = &coordTrue[1][2][0];
        coord[2][0] = &coordTrue[2][0][0]; coord[2][1] = &coordTrue[2][1][0];
        coord[2][2] = &coordTrue[2][2][0];
    }

    PetscFunctionReturn(0);
}


// create string containing information of this mesh
PetscErrorCode CartesianMesh::createInfoString()
{
    PetscFunctionBeginUser;

    PetscErrorCode          ierr;

    std::stringstream       ss;


    if (mpiRank == 0)
    {
        ss << std::string(80, '=') << std::endl;
        ss << "Cartesian Staggered Grid:" << std::endl;
        ss << std::string(80, '=') << std::endl;


        ss << "\tDimension: " << dim << std::endl;
        ss << std::endl;

        ss << "\tDomain Range: " 
            << "[" << min[0] << ", " << max[0] << "]; "
            << "[" << min[1] << ", " << max[1] << "]";
        if (dim == 3) ss << "; [" << min[2] << ", " << max[2] << "]";
        ss << std::endl;
        ss << std::endl;

        ss << "\tNumber of Cells (Nx x Ny" << ((dim==2)? "" : " x Nz") << "):\n";
        ss << "\t\tPressure Grid: ";
        ss << n[3][0];
        for(int dir=1; dir<dim; ++dir) ss << " x " << n[3][dir];
        ss << std::endl;

        for(int comp=0; comp<dim; ++comp)
        {
            ss << "\t\t" << fd2str[Field(comp)] 
                << "-Velocity Grid: ";
            ss << n[comp][0];
            for(int dir=1; dir<dim; ++dir) ss << " x " << n[comp][dir];
            ss << std::endl;
        }
        ss << std::endl;

        ss << "\tGrid Decomposition:" << std::endl;
        ss << "\t\tMPI Processes: " 
            << nProc[0] << " x " << nProc[1] << " x " << nProc[2] << std::endl;
    }

    ierr = addLocalInfoString(ss); CHKERRQ(ierr);
    ss << std::endl;

    info = ss.str();

    PetscFunctionReturn(0);
}


// create local information of mesh on this process
PetscErrorCode CartesianMesh::addLocalInfoString(std::stringstream &ss)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    std::string     pre("\t\t");

    ss << pre << "Rank " << mpiRank << ": " << std::endl;
    ss << pre << "\tPressure Grid: ";
    ss << "[" << bg[3][0] << ", " << ed[3][0] << "), ";
    ss << "[" << bg[3][1] << ", " << ed[3][1] << "), ";
    ss << "[" << bg[3][2] << ", " << ed[3][2] << ") ";
    ss << std::endl;

    for(int comp=0; comp<dim; ++comp)
    {
        ss << pre << "\t" << fd2str[Field(comp)] << "-Velocity Grid: ";
        ss << "[" << bg[comp][0] << ", " << ed[comp][0] << "), ";
        ss << "[" << bg[comp][1] << ", " << ed[comp][1] << "), ";
        ss << "[" << bg[comp][2] << ", " << ed[comp][2] << ") ";
        ss << std::endl;
    }

    PetscFunctionReturn(0);
}


// initialize PETSc DMDA objects
PetscErrorCode CartesianMesh::initDMDA(const BoolVec2D &periodic)
{
    PetscFunctionBeginUser;

    PetscErrorCode  ierr;

    ierr = createPressureDMDA(periodic); CHKERRQ(ierr);
    ierr = createVelocityPack(periodic); CHKERRQ(ierr);

    // gather numbers of local velocity component
    for(PetscInt f=0; f<dim; ++f)
    {
        // temporarily use UNLocal
        UNLocal = m[f][0] * m[f][1] * m[f][2];

        // get number of local velocity component from all processes
        ierr = MPI_Barrier(comm); CHKERRQ(ierr);
        ierr = MPI_Allgather(&UNLocal, 1, MPIU_INT,
                UNLocalAllProcs[f].data(), 1, MPIU_INT, comm); CHKERRQ(ierr);
    }

    // total number of local velocity points
    UNLocal = UNLocalAllProcs[0][mpiRank] + 
        UNLocalAllProcs[1][mpiRank] + UNLocalAllProcs[2][mpiRank];

    // offset on each process of each velocity component
    for(PetscInt f=0; f<dim; ++f)
    {
        for(PetscMPIInt r=mpiSize-1; r>0; --r)
            offsetsAllProcs[f][r] = UNLocalAllProcs[f][r-1];

        for(PetscMPIInt r=1; r<mpiSize; ++r)
            offsetsAllProcs[f][r] += offsetsAllProcs[f][r-1];
    }

    // calculate total number of all local velocity points on all processes
    for(PetscMPIInt r=0; r<mpiSize; ++r) UPackNLocalAllProcs[r] = 
        UNLocalAllProcs[0][r] + UNLocalAllProcs[1][r] + UNLocalAllProcs[2][r];

    // offsets of all velocity points on each process in packed form
    for(PetscMPIInt r=mpiSize-1; r>0; --r)
        offsetsPackAllProcs[r] = UPackNLocalAllProcs[r-1];

    for(PetscMPIInt r=1; r<mpiSize; ++r)
        offsetsPackAllProcs[r] += offsetsPackAllProcs[r-1];


    // total number of local pressure points
    pNLocal = m[3][0] * m[3][1] * m[3][2];

    PetscFunctionReturn(0);
}


// create a single PETSc DMDA object
PetscErrorCode CartesianMesh::createSingleDMDA(
        const PetscInt &i, const BoolVec2D &periodic)
{
    using namespace std;

    PetscFunctionBeginUser;

    PetscErrorCode  ierr;

    DMDALocalInfo   lclInfo;
    

    switch (dim)
    {
        case 2:
            ierr = DMDACreate2d(comm,
                    (periodic[0][0])? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED,
                    (periodic[1][1])? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED,
                    DMDA_STENCIL_STAR, 
                    n[i][0], n[i][1], nProc[0], nProc[1],
                    1, 1, nullptr, nullptr, &da[i]); CHKERRQ(ierr);
            break;
        case 3:
            ierr = DMDACreate3d(comm, 
                    (periodic[0][0])? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED,
                    (periodic[1][1])? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED,
                    (periodic[2][2])? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED,
                    DMDA_STENCIL_STAR, n[i][0], n[i][1], n[i][2], 
                    nProc[0], nProc[1], nProc[2], 
                    1, 1, nullptr, nullptr, nullptr, &da[i]); CHKERRQ(ierr);
            break;
    }

    ierr = DMDAGetLocalInfo(da[i], &lclInfo); CHKERRQ(ierr);

    bg[i][0] = lclInfo.xs; bg[i][1] = lclInfo.ys; bg[i][2] = lclInfo.zs;

    m[i][0] = lclInfo.xm; m[i][1] = lclInfo.ym; m[i][2] = lclInfo.zm;

    for(unsigned int comp=0; comp<3; ++comp) 
        ed[i][comp] = bg[i][comp] + m[i][comp];

    PetscFunctionReturn(0);
}


// create a PETSc DMDA for pressure field
PetscErrorCode CartesianMesh::createPressureDMDA(const BoolVec2D &periodic)
{
    using namespace std;

    PetscFunctionBeginUser;

    PetscErrorCode  ierr;

    ierr = createSingleDMDA(3, periodic); CHKERRQ(ierr);
    
    ierr = DMDAGetAO(da[3], &ao[3]); CHKERRQ(ierr);

    ierr = DMDAGetInfo(da[3], nullptr, nullptr, nullptr, nullptr,
           &nProc[0], &nProc[1], &nProc[2], 
           nullptr, nullptr, nullptr, nullptr, nullptr, nullptr); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


// create a PETSc DMComposite object for velocity fields
PetscErrorCode CartesianMesh::createVelocityPack(const BoolVec2D &periodic)
{
    using namespace std;

    PetscFunctionBeginUser;

    PetscErrorCode  ierr;

    ierr = DMCompositeCreate(comm, &UPack); CHKERRQ(ierr);

    for(int i=0; i<dim; ++i)
    {
        ierr = createSingleDMDA(i, periodic); CHKERRQ(ierr);
        ierr = DMDAGetAO(da[i], &ao[i]); CHKERRQ(ierr);
        ierr = DMCompositeAddDM(UPack, da[i]); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}


// get natural index through MatStencil
PetscErrorCode CartesianMesh::getNaturalIndex(
        const PetscInt &f, const MatStencil &s, PetscInt &idx) const
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    ierr = getNaturalIndex(f, s.i, s.j, s.k, idx); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


// get natural index through i, j, k indices
PetscErrorCode CartesianMesh::getNaturalIndex(const PetscInt &f, 
        const PetscInt &i, const PetscInt &j, const PetscInt &k, 
        PetscInt &idx) const
{
    PetscFunctionBeginUser;

    // some overhead due to implicit if-condition, but this is also how PETSc
    // does in their source code for similar functions
    idx = i + j * n[f][0] + ((dim == 3) ? (k * n[f][1] * n[f][0]) : 0);

    PetscFunctionReturn(0);
}


// get local index through MatStencil
PetscErrorCode CartesianMesh::getLocalIndex(
        const PetscInt &f, const MatStencil &s, PetscInt &idx) const
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    ierr = DMDAConvertToCell(da[f], s, &idx); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


// get local index through i, j, k indices
PetscErrorCode CartesianMesh::getLocalIndex(const PetscInt &f, 
        const PetscInt &i, const PetscInt &j, const PetscInt &k, 
        PetscInt &idx) const
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    ierr = getLocalIndex(f, {k, j, i, 0}, idx); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


// get global index through MatStencil
PetscErrorCode CartesianMesh::getGlobalIndex(
        const PetscInt &f, const MatStencil &s, PetscInt &idx) const
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    ierr = getNaturalIndex(f, s, idx); CHKERRQ(ierr);
    ierr = AOApplicationToPetsc(ao[f], 1, &idx); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


// get global index through i, j, k indices
PetscErrorCode CartesianMesh::getGlobalIndex(const PetscInt &f, 
        const PetscInt &i, const PetscInt &j, const PetscInt &k, 
        PetscInt &idx) const
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    ierr = getGlobalIndex(f, {k, j, i, 0}, idx); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


// get index in DMComposite through MatStencil
PetscErrorCode CartesianMesh::getPackedGlobalIndex(
        const PetscInt &f, const MatStencil &s, PetscInt &idx) const
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    PetscMPIInt         p;

    PetscInt            unPackedIdx;

    // get the global index of the target in sub-DM (un-packed DM)
    ierr = getGlobalIndex(f, s, unPackedIdx); CHKERRQ(ierr);
    
    if (f == 3) // pressure isn't in packed DM, so packed-ID = global-ID.
    {
        idx = unPackedIdx;
        PetscFunctionReturn(0);
    }

    // find which process owns the target
    p = std::upper_bound(offsetsAllProcs[f].begin(), offsetsAllProcs[f].end(), 
            unPackedIdx) - offsetsAllProcs[f].begin() - 1;

    // the beginngin index of process p in packed DM
    idx = offsetsPackAllProcs[p]; 

    // offset due to previous DMs
    for(PetscInt i=0; i<f; ++i) idx += UNLocalAllProcs[i][p]; // offset from previous DMs

    // the packed index will be this
    idx += (unPackedIdx - offsetsAllProcs[f][p]);

    PetscFunctionReturn(0);
}


// get indec in DMComposite through i, j, k indices
PetscErrorCode CartesianMesh::getPackedGlobalIndex(const PetscInt &f, 
        const PetscInt &i, const PetscInt &j, const PetscInt &k, 
        PetscInt &idx) const
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    ierr = getPackedGlobalIndex(f, {k, j, i, 0}, idx); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

} // end of namespace mesh
} // end of namespace petibm
