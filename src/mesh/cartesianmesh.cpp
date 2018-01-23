/**
 * \file cartesianmesh.cpp
 * \brief Implementations of mesh::CartesianMesh.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
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
# include <petibm/io.h>


namespace petibm
{
namespace mesh
{

using namespace type;


// implementation of CartesianMesh::CastesianMesh
CartesianMesh::CartesianMesh(const MPI_Comm &world, const YAML::Node &node)
{
    init(world, node);
} // CartesianMesh


// implementation of CartesianMesh::~CastesianMesh
CartesianMesh::~CartesianMesh()
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;
    PetscBool finalized;

    ierr = PetscFinalized(&finalized); CHKERRV(ierr);
    if (finalized) return;

    std::vector<AO>().swap(ao); // DO NOT DESTROY UNDERLYING AOs!!
} // ~CartesianMesh


// implementation of CartesianMesh::destroy
PetscErrorCode CartesianMesh::destroy()
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    type::RealVec3D().swap(dLTrue);
    type::RealVec3D().swap(coordTrue);
    std::vector<AO>().swap(ao); // DO NOT DESTROY UNDERLYING AOs!!
    type::IntVec1D().swap(UPackNLocalAllProcs);
    type::IntVec2D().swap(offsetsAllProcs);
    type::IntVec1D().swap(offsetsPackAllProcs);

    ierr = MeshBase::destroy(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // destroy


// implementation of CartesianMesh::init
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
    // with carefully chosen default values, we can take care of 2D case with 
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
    ierr = parser::parseBCs(node, bcTypes, bcValues); CHKERRQ(ierr);
    ierr = misc::checkPeriodicBC(bcTypes, periodic); CHKERRQ(ierr);
    
    // create raw grid information
    ierr = createPressureMesh(); CHKERRQ(ierr);
    ierr = createVertexMesh(); CHKERRQ(ierr);
    ierr = createVelocityMesh(); CHKERRQ(ierr);
    ierr = MPI_Barrier(comm); CHKERRQ(ierr);

    // create PETSc DMs
    ierr = initDMDA(); CHKERRQ(ierr);

    // setup the format of the file that the `write` function will use
    //ierr = setOutputFormat(type); CHKERRQ(ierr);

    // create a std::string that can be used in `printInfo` and output stream
    ierr = createInfoString(); CHKERRQ(ierr);

    // all processes should be synchronized
    ierr = MPI_Barrier(comm); CHKERRQ(ierr);


    PetscFunctionReturn(0);
} // init


// implementation of CartesianMesh::createPressureMesh
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

    // in 2D simulation, we still use the variable dL in z-direction
    if (dim == 2)
    {
        dL[3][2] = &dLTrue[3][2][0];
        coord[3][2] = &coordTrue[3][2][0];
    }

    PetscFunctionReturn(0);
} // createPressureMesh


// implementation of CartesianMesh::createVertexMesh
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

        // shift according to real coord[4][i][0], which is min
        std::for_each(coordTrue[4][i].begin(), coordTrue[4][i].end(),
                [i, this] (double &a) -> void { a += this->min[i]; });

        // no ghost points in vertex mesh, so coord[4][i][0]=coordTrue[4][i][0]
        coord[4][i] = &coordTrue[4][i][0];
    }

    // in 2D simulation, we still use the variable coord in z-direction
    if (dim == 2) coord[4][2] = &coordTrue[4][2][0];

    PetscFunctionReturn(0);
} // createVertexMesh


// implementation of CartesianMesh::createCelocityMesh
PetscErrorCode CartesianMesh::createVelocityMesh()
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

                // coordinate and cell size will match that of pressure point,
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
} // createVelocityMesh


// implementation of CartesianMesh::createInfoString
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
} // createInfoString


// implementation of CartesianMesh::addLocalInfoString
PetscErrorCode CartesianMesh::addLocalInfoString(std::stringstream &ss)
{
    PetscFunctionBeginUser;

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
} // addLocalInfoString


// implementation of CartesianMesh::initDMDA
PetscErrorCode CartesianMesh::initDMDA()
{
    PetscFunctionBeginUser;

    PetscErrorCode  ierr;

    ierr = createPressureDMDA(); CHKERRQ(ierr);
    ierr = createVelocityPack(); CHKERRQ(ierr);

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
} // initDMDA


// implementation of CartesianMesh::createSingleDMDA
PetscErrorCode CartesianMesh::createSingleDMDA(const PetscInt &i)
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
                    DMDA_STENCIL_BOX, 
                    n[i][0], n[i][1], nProc[0], nProc[1],
                    1, 1, nullptr, nullptr, &da[i]); CHKERRQ(ierr);
            break;
        case 3:
            ierr = DMDACreate3d(comm, 
                    (periodic[0][0])? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED,
                    (periodic[1][1])? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED,
                    (periodic[2][2])? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_GHOSTED,
                    DMDA_STENCIL_BOX, n[i][0], n[i][1], n[i][2], 
                    nProc[0], nProc[1], nProc[2], 
                    1, 1, nullptr, nullptr, nullptr, &da[i]); CHKERRQ(ierr);
            break;
    }
    
    ierr = DMSetUp(da[i]); CHKERRQ(ierr);

    ierr = DMDAGetLocalInfo(da[i], &lclInfo); CHKERRQ(ierr);

    bg[i][0] = lclInfo.xs; bg[i][1] = lclInfo.ys; bg[i][2] = lclInfo.zs;

    m[i][0] = lclInfo.xm; m[i][1] = lclInfo.ym; m[i][2] = lclInfo.zm;

    for(unsigned int comp=0; comp<3; ++comp) 
        ed[i][comp] = bg[i][comp] + m[i][comp];

    PetscFunctionReturn(0);
} // createSingleDMDA


// implementation of CartesianMesh::createPressureDMDA
PetscErrorCode CartesianMesh::createPressureDMDA()
{
    using namespace std;

    PetscFunctionBeginUser;

    PetscErrorCode  ierr;

    ierr = createSingleDMDA(3); CHKERRQ(ierr);
    
    ierr = DMDAGetAO(da[3], &ao[3]); CHKERRQ(ierr);

    ierr = DMDAGetInfo(da[3], nullptr, nullptr, nullptr, nullptr,
           &nProc[0], &nProc[1], &nProc[2], 
           nullptr, nullptr, nullptr, nullptr, nullptr, nullptr); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // createPressureDMDA


// implementation of CartesianMesh::createVelocityPack
PetscErrorCode CartesianMesh::createVelocityPack()
{
    using namespace std;

    PetscFunctionBeginUser;

    PetscErrorCode  ierr;

    ierr = DMCompositeCreate(comm, &UPack); CHKERRQ(ierr);

    for(int i=0; i<dim; ++i)
    {
        ierr = createSingleDMDA(i); CHKERRQ(ierr);
        ierr = DMDAGetAO(da[i], &ao[i]); CHKERRQ(ierr);
        ierr = DMCompositeAddDM(UPack, da[i]); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
} // createVelocityPack


// implementation of CartesianMesh::getNaturalIndex
PetscErrorCode CartesianMesh::getNaturalIndex(
        const PetscInt &f, const MatStencil &s, PetscInt &idx) const
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    ierr = getNaturalIndex(f, s.i, s.j, s.k, idx); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // getNaturalIndex


// implementation of CartesianMesh::getNaturalIndex
PetscErrorCode CartesianMesh::getNaturalIndex(const PetscInt &f, 
        const PetscInt &i, const PetscInt &j, const PetscInt &k, 
        PetscInt &idx) const
{
    PetscFunctionBeginUser;

# ifndef NDEBUG
    if ((i < -1) || (i > n[f][0]) ||
        (j < -1) || (j > n[f][1]) ||
        (k < -1) || (k > n[f][2]))
        SETERRQ4(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
           "Stencil (%d, %d, %d) of field %d is out of domain.\n",
           i, j, k, f);
# endif
    
    if (i == -1)
    {
        if (periodic[0][0])
            idx = (n[f][0] - 1) + j * n[f][0] + ((dim == 3) ? (k * n[f][1] * n[f][0]) : 0);
        else
            idx = -1;
        
        PetscFunctionReturn(0);
    }
    
    if (i == n[f][0])
    {
        if (periodic[0][0])
            idx = j * n[f][0] + ((dim == 3) ? (k * n[f][1] * n[f][0]) : 0);
        else
            idx = -1;
        
        PetscFunctionReturn(0);
    }
    
    if (j == -1)
    {
        if (periodic[0][1])
            idx = i + (n[f][1]-1) * n[f][0] + ((dim == 3) ? (k * n[f][1] * n[f][0]) : 0);
        else
            idx = -1;
        
        PetscFunctionReturn(0);
    }
    
    if (j == n[f][1])
    {
        if (periodic[0][1])
            idx = i + ((dim == 3) ? (k * n[f][1] * n[f][0]) : 0);
        else
            idx = -1;
        
        PetscFunctionReturn(0);
    }
    
    if (k == -1)
    {
        if (periodic[0][2])
            idx = i + j * n[f][0] + ((dim == 3) ? ((n[f][2] - 1) * n[f][1] * n[f][0]) : 0);
        else
            idx = -1;
        
        PetscFunctionReturn(0);
    }
    
    if (k == n[f][2])
    {
        if (periodic[0][2])
            idx = i + j * n[f][0];
        else
            idx = -1;
        
        PetscFunctionReturn(0);
    }

    // some overhead due to implicit if-condition, but this is also how PETSc
    // does in their source code for similar functions
    idx = i + j * n[f][0] + ((dim == 3) ? (k * n[f][1] * n[f][0]) : 0);

    PetscFunctionReturn(0);
} // getNaturalIndex


// implementation of CartesianMesh::getLocalIndex
PetscErrorCode CartesianMesh::getLocalIndex(
        const PetscInt &f, const MatStencil &s, PetscInt &idx) const
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    ierr = DMDAConvertToCell(da[f], s, &idx); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // getLocalIndex


// implementation of CartesianMesh::getLocalIndex
PetscErrorCode CartesianMesh::getLocalIndex(const PetscInt &f, 
        const PetscInt &i, const PetscInt &j, const PetscInt &k, 
        PetscInt &idx) const
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    ierr = getLocalIndex(f, {k, j, i, 0}, idx); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // getLocalIndex


// implementation of CartesianMesh::getGlobalIndex
PetscErrorCode CartesianMesh::getGlobalIndex(
        const PetscInt &f, const MatStencil &s, PetscInt &idx) const
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    ierr = getNaturalIndex(f, s, idx); CHKERRQ(ierr);
    ierr = AOApplicationToPetsc(ao[f], 1, &idx); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // getGlobalIndex


// implementation of CartesianMesh::getGlobalIndex
PetscErrorCode CartesianMesh::getGlobalIndex(const PetscInt &f, 
        const PetscInt &i, const PetscInt &j, const PetscInt &k, 
        PetscInt &idx) const
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    ierr = getGlobalIndex(f, {k, j, i, 0}, idx); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // getGlobalIndex


// implementation of CartesianMesh::getPackedGlobalIndex
PetscErrorCode CartesianMesh::getPackedGlobalIndex(
        const PetscInt &f, const MatStencil &s, PetscInt &idx) const
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    PetscMPIInt         p;

    PetscInt            unPackedIdx;

    // get the global index of the target in sub-DM (un-packed DM)
    ierr = getGlobalIndex(f, s, unPackedIdx); CHKERRQ(ierr);
    
    // 1. pressure isn't in packed DM, so packed-ID = global-ID.
    // 2. if unPackedIdx is -1, then packed ID is -1, too.
    if ((f == 3) || (unPackedIdx == -1))
    {
        idx = unPackedIdx;
        PetscFunctionReturn(0);
    }

    // find which process owns the target
    p = std::upper_bound(offsetsAllProcs[f].begin(), offsetsAllProcs[f].end(), 
            unPackedIdx) - offsetsAllProcs[f].begin() - 1;

    // the beginning index of process p in packed DM
    idx = offsetsPackAllProcs[p]; 

    // offset due to previous DMs
    for(PetscInt i=0; i<f; ++i) idx += UNLocalAllProcs[i][p]; // offset from previous DMs

    // the packed index will be this
    idx += (unPackedIdx - offsetsAllProcs[f][p]);

    PetscFunctionReturn(0);
} // getPackedGlobalIndex


// implementation of CartesianMesh::getPackedGlobalIndex
PetscErrorCode CartesianMesh::getPackedGlobalIndex(const PetscInt &f, 
        const PetscInt &i, const PetscInt &j, const PetscInt &k, 
        PetscInt &idx) const
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    ierr = getPackedGlobalIndex(f, {k, j, i, 0}, idx); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // getPackedGlobalIndex


// implementation of CartesianMesh::write
PetscErrorCode CartesianMesh::write(const std::string &filePath) const
{
    PetscFunctionBeginUser;
    
    PetscErrorCode      ierr;
    
    if (mpiRank == 0)
    {
        PetscFileMode   mode = FILE_MODE_WRITE;
        std::vector<std::string>    names{"x", "y", "z"};
        
        for(unsigned int f=0; f<5; ++f)
        {
            std::string group = type::fd2str[type::Field(f)];
            
            ierr = io::writeHDF5Vecs(PETSC_COMM_SELF, filePath, group,
                    names, n[f], coord[f], mode); CHKERRQ(ierr);
            
            mode = FILE_MODE_APPEND;
        }
    }
    
    ierr = MPI_Barrier(comm); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
} // write

} // end of namespace mesh
} // end of namespace petibm
