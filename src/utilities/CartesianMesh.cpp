/***************************************************************************//**
 * \file CartesianMesh.cpp
 * \author Anush Krishnan (anush@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Implementation of the methods of the class `CartesianMesh`.
 */


// here goes C++ STL
# include <iostream>
# include <fstream>
# include <sstream>
# include <iomanip>
# include <algorithm>
# include <numeric>

// here goes PETSc headers
# include <petscvec.h>
# include <petscviewerhdf5.h>

// here goes HDF5
# include "hdf5.h"

// here goes headers from our PetIBM
# include "CartesianMesh.h"
# include "parser.h"


using namespace types;


//  TODO: move CMD parsing to other part in PetIBM
CartesianMesh::CartesianMesh() = default;


CartesianMesh::CartesianMesh(
        const MPI_Comm &world, const YAML::Node &node, 
        BCInfoHolder &bc, const OutputType &type)
{
    init(world, node, bc, type);
}


CartesianMesh::~CartesianMesh() = default;


PetscErrorCode CartesianMesh::init(
        const MPI_Comm &world, const YAML::Node &meshNode, 
        BCInfoHolder &bc, const OutputType &type)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    // store the address of the communicator
    // note: this is a bad practice; shared_ptr is not for stack variables!!
    comm = std::shared_ptr<const MPI_Comm>(&world, [](const MPI_Comm*){});

    // set rank and size
    ierr = MPI_Comm_size(*comm, &mpiSize); CHKERRQ(ierr);
    ierr = MPI_Comm_rank(*comm, &mpiRank); CHKERRQ(ierr);


    // set the default sizes and values for members
    // The default values are chosen for 2D cases. We want to eliminate the use
    // of template and also `if` conditions that determine the dimension. So
    // with carefully choosed default values, we can take care of 2D case with 
    // 3D code in many places.
    min = types::RealVec1D(3, 0.0);
    max = types::RealVec1D(3, 1.0);
    n = types::IntVec2D(5, types::IntVec1D(3, 1));
    coord = types::RealVec3D(5, types::RealVec2D(3, types::RealVec1D(1, 0.0)));
    dLTrue = types::RealVec3D(5, types::RealVec2D(3, types::RealVec1D(1, 1.0)));
    dL = types::DeltaLVec(5, std::vector<PetscReal*>(3, nullptr));
    da = std::vector<DM>(5, PETSC_NULL);
    nProc = types::IntVec1D(3, PETSC_DECIDE);
    bg = types::IntVec2D(5, types::IntVec1D(3, 0));
    ed = types::IntVec2D(5, types::IntVec1D(3, 1));
    m = types::IntVec2D(5, types::IntVec1D(3, 0));


    // store the address of bcInfo
    // note: this is a bad practice; shared_ptr is not for stack variables!!
    bcInfo = std::shared_ptr<BCInfoHolder>(&bc, [](BCInfoHolder*){});


    // index 3 represent pressure mesh; min & max always represent pressure mesh
    ierr = parser::parseMesh(meshNode, dim, min, max, n[3], dLTrue[3]); CHKERRQ(ierr);


    // create raw grid information
    ierr = createPressureMesh(); CHKERRQ(ierr);
    ierr = createVertexMesh(); CHKERRQ(ierr);
    ierr = createVelocityMesh(); CHKERRQ(ierr);
    ierr = MPI_Barrier(*comm); CHKERRQ(ierr);

    // create PETSc DMs
    ierr = initDMDA(); CHKERRQ(ierr);

    // setup the format of the file that the `write` function will use
    ierr = setOutputFormat(type); CHKERRQ(ierr);

    // create a std::string that can be used in `printInfo` and output stream
    ierr = createInfoString(); CHKERRQ(ierr);

    // all processes should be syncrinized
    ierr = MPI_Barrier(*comm); CHKERRQ(ierr);


    PetscFunctionReturn(0);
}


PetscErrorCode CartesianMesh::setOutputFormat(const OutputType &type)
{
    PetscFunctionBeginUser;

    using namespace std::placeholders;

    switch (type)
    {
        case OutputType::Binary:
            write = std::bind(&CartesianMesh::writeBinary, this, _1, _2);
            break;
        case OutputType::VTK:
            write = std::bind(&CartesianMesh::writeVTK, this, _1, _2);
            break;
        case OutputType::HDF5:
            write = std::bind(&CartesianMesh::writeHDF5, this, _1, _2);
            break;
    }

    PetscFunctionReturn(0);
}


PetscErrorCode CartesianMesh::createPressureMesh()
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    // loop through all axes to get the coordinates of pressure points
    for(unsigned int i=0; i<dim; ++i)
    {
        coord[3][i].resize(n[3][i]);

        // the following three lines calculate the coord of pressure points
        std::partial_sum(
                dLTrue[3][i].begin(), dLTrue[3][i].end(), coord[3][i].begin());

        auto f = [i, this] (double &a, double &b) -> double 
                    { return a + this->min[i] - 0.5 * b; };

        std::transform(coord[3][i].begin(), coord[3][i].end(), 
                dLTrue[3][i].begin(), coord[3][i].begin(), f);

        // no ghost points in pressure grid, so dL[0] = dLTrue[0]
        dL[3][i] = &dLTrue[3][i][0];
    }

    // in 2D simulation, we still use the varaible dL in z-direction
    if (dim == 2) dL[3][2] = &dLTrue[3][2][0];

    PetscFunctionReturn(0);
}


PetscErrorCode CartesianMesh::createVertexMesh()
{
    PetscFunctionBeginUser;

    // loop through all axes to get the coordinates of mesh vertexes
    // note: index 3 means pressure mesh; 4 means vertexes
    for(unsigned int i=0; i<dim; ++i)
    {
        // number of vertexes is one more than pressure cells
        n[4][i] = n[3][i] + 1;
        coord[4][i].resize(n[4][i]);

        // create coordinates of vertexes
        std::partial_sum(dLTrue[3][i].begin(), dLTrue[3][i].end(), 
                coord[4][i].begin()+1); // here we assume coord[4][0] = 0.0

        std::for_each(coord[4][i].begin(), coord[4][i].end(),
                [i, this] (double &a) -> void { a += this->min[i]; });
    }

    PetscFunctionReturn(0);
}


PetscErrorCode CartesianMesh::createVelocityMesh()
{
    PetscFunctionBeginUser;

    // loop through all velocity component
    // note: index 0, 1, 2 represent u, v, w respectively
    for(unsigned int comp=0; comp<dim; comp++)
    {
        // loop through each direction, i.e., x, y, z directions
        for(unsigned int dir=0; dir<dim; ++dir)
        {
            // when the direction corresponding to velocity component
            if (dir == comp)
            {
                // there will no velocity point on boundary (it's a ghost point)
                n[comp][dir] = n[3][dir] - 1;

                // coordinates will match the interior part of the vertex grid
                coord[comp][dir] = RealVec1D(
                        coord[4][dir].begin()+1, coord[4][dir].end()-1);

                // we store the dL of ghost points, so the size is n+2
                dLTrue[comp][dir].resize(n[comp][dir]+2);

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
                if ((*bcInfo)[BCLoc(dir*2)][Field(comp)].type == BCType::PERIODIC)
                {
                    // add 1 to the number of valid grid point
                    n[comp][dir] += 1;

                    // add the coordinate to the extra point
                    coord[comp][dir].push_back(coord[4][dir].back());


                    // the space for larger ghost point in dLTrue is now used to 
                    // store the grid point on periodic BC
                    dLTrue[comp][dir].back() = 
                        f(dLTrue[3][dir][0], dLTrue[3][dir].back());

                    // the smaller ghost point will be the same as the last
                    // valid point
                    dLTrue[comp][dir][0] = dLTrue[comp][dir].back();

                    // we have to create one more space for the larger ghost
                    // point, and its dL is the same as the 1st valid point
                    dLTrue[comp][dir].push_back(dLTrue[comp][dir][1]);
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

                // coordinate and cel size will match that of pressure point
                coord[comp][dir] = coord[3][dir];

                // we store the dL of ghost points, so the size is n+2
                dLTrue[comp][dir].resize(n[comp][dir]+2);

                // the dL of valid points is the same as pressure grid
                std::copy(dLTrue[3][dir].begin(), dLTrue[3][dir].end(), 
                        dLTrue[comp][dir].begin()+1);

                // get the dL for the smaller and larger ghost points
                dLTrue[comp][dir].front() = dLTrue[3][dir].front();
                dLTrue[comp][dir].back() = dLTrue[3][dir].back();
            }


            // variable dL will point to the 1st valid grid point, so we can access
            // the 1st ghost point with the index -1, i.e., dL[comp][dir][-1]
            dL[comp][dir] = &dLTrue[comp][dir][1];
        }
    }


    // for 2D simulation, we still use dL in z-direction in our code
    if (dim == 2)
    {
        dL[0][2] = &dLTrue[0][2][0];
        dL[1][2] = &dLTrue[1][2][0];
        dL[2][0] = &dLTrue[2][0][0];
        dL[2][1] = &dLTrue[2][1][0];
        dL[2][2] = &dLTrue[2][2][0];
    }

    PetscFunctionReturn(0);
}


PetscErrorCode CartesianMesh::createInfoString()
{
    PetscFunctionBeginUser;

    PetscErrorCode          ierr;

    std::stringstream       ss;


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
    for(unsigned int dir=1; dir<dim; ++dir) ss << " x " << n[3][dir];
    ss << std::endl;

    for(unsigned int comp=0; comp<dim; ++comp)
    {
        ss << "\t\t" << fd2str[Field(comp)] 
            << "-Velocity Grid: ";
        ss << n[comp][0];
        for(unsigned int dir=1; dir<dim; ++dir) ss << " x " << n[comp][dir];
        ss << std::endl;
    }
    ss << std::endl;

    ss << "\tGrid Decomposition:" << std::endl;
    ss << "\t\tMPI Processes: " 
        << nProc[0] << " x " << nProc[1] << " x " << nProc[2] << std::endl;

    ierr = addLocalInfoString(ss); CHKERRQ(ierr);
    ss << std::endl;

    info = ss.str();

    PetscFunctionReturn(0);
}


PetscErrorCode CartesianMesh::addLocalInfoString(std::stringstream &ss)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    std::string     pre("\t\t");

    IntVec2D        start = IntVec2D(5, IntVec1D(mpiSize * 3)),
                    end = IntVec2D(5, IntVec1D(mpiSize * 3));

    for(unsigned int i=0; i<dim; ++i)
    {
        ierr = MPI_Allgather(bg[i].data(), 3, MPIU_INT, 
                start[i].data(), 3, MPIU_INT, *comm); CHKERRQ(ierr);
        ierr = MPI_Allgather(ed[i].data(), 3, MPIU_INT, 
                end[i].data(), 3, MPIU_INT, *comm); CHKERRQ(ierr);
    }

    ierr = MPI_Allgather(bg[3].data(), 3, MPIU_INT, 
            start[3].data(), 3, MPIU_INT, *comm); CHKERRQ(ierr);
    ierr = MPI_Allgather(ed[3].data(), 3, MPIU_INT, 
            end[3].data(), 3, MPIU_INT, *comm); CHKERRQ(ierr);


    for(PetscMPIInt i=0; i<mpiSize; ++i)
    {
        ss << pre << "Rank " << i << ": " << std::endl;
        ss << pre << "\tPressure Grid: ";
        ss << "[" << start[3][i*3] << ", " << end[3][i*3] << "), ";
        ss << "[" << start[3][i*3+1] << ", " << end[3][i*3+1] << "), ";
        ss << "[" << start[3][i*3+2] << ", " << end[3][i*3+2] << ") ";
        ss << std::endl;

        for(unsigned int comp=0; comp<dim; ++comp)
        {
            ss << pre << "\t" << fd2str[Field(comp)] << "-Velocity Grid: ";
            ss << "[" << start[comp][i*3] << ", " << end[comp][i*3] << "), ";
            ss << "[" << start[comp][i*3+1] << ", " << end[comp][i*3+1] << "), ";
            ss << "[" << start[comp][i*3+2] << ", " << end[comp][i*3+2] << ") ";
            ss << std::endl;
        }
    }

    ierr = MPI_Barrier(*comm); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


PetscErrorCode CartesianMesh::printInfo() const
{
    PetscFunctionBeginUser;
    
    PetscErrorCode ierr;
    ierr = PetscPrintf(*comm, info.c_str()); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // printInfo


PetscErrorCode CartesianMesh::writeBinary(
        const std::string &dir, const std::string &file)
{
    PetscFunctionBeginUser;

    PetscErrorCode ierr;

    if (mpiRank == 0)
    {
        std::ofstream   streamFile(dir + "/" + file + ".txt");

        streamFile << n[3][0] << "\t" << n[3][1] << "\t"
            << ((dim==3)? std::to_string(n[3][2]) : "") << std::endl;

        for(auto it: coord[4][0]) 
            streamFile << std::setprecision(16) << it << std::endl;

        for(auto it: coord[4][1]) 
            streamFile << std::setprecision(16) << it << std::endl;

        if (dim == 3)
            for(auto it: coord[4][2]) 
                streamFile << std::setprecision(16) << it << std::endl;

        streamFile.close();
    }

    ierr = MPI_Barrier(*comm); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // write


PetscErrorCode CartesianMesh::writeVTK(
        const std::string &dir, const std::string &file)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    // only the master process will write file
    if (mpiRank == 0)
    {
        std::ofstream   fs(dir + "/" + file + ".vtk");

        // vtk version
        fs << "# vtk DataFile Version 3.0" << std::endl;

        // one-line comment about this grid file
        fs << "Cartesian Grid" << std::endl;

        // so far, we use ASCII for this vtk file
        fs << "ASCII" << std::endl;

        // use RECTILINEAR_GRID to save the file size
        fs << "DATASET RECTILINEAR_GRID" << std::endl;

        // dimension of the vertexes (not the pressure cells!)
        fs << "DIMENSIONS ";
        for(auto it: n[4]) fs << it << " ";
        fs << std::endl;

        // x coordinates of the vertexes
        fs << "X_COORDINATES " << n[4][0] << " double" << std::endl;
        for(auto it: coord[4][0]) fs << std::setprecision(16) << it << " ";
        fs << std::endl;

        // y coordinates of the vertexes
        fs << "Y_COORDINATES " << n[4][1] << " double" << std::endl;
        for(auto it: coord[4][1]) fs << std::setprecision(16) << it << " ";
        fs << std::endl;

        // z coordinates of the vertexes
        fs << "Z_COORDINATES " << n[4][2] << " double" << std::endl;
        for(auto it: coord[4][2]) fs << std::setprecision(16) << it << " ";
        fs << std::endl;

        // close the file
        fs.close();
    }

    // all processes must to wait the master process
    ierr = MPI_Barrier(*comm); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


PetscErrorCode CartesianMesh::writeHDF5(
        const std::string &dir, const std::string &file)
{
    PetscFunctionBeginUser;

# ifndef PETSC_HAVE_HDF5
    SETERRQ(*comm, 56, "Seems the PETSc was not compiled with HDF5.");
# else

    PetscErrorCode  ierr;

    if (mpiRank == 0)
    {
        std::string     extFile = dir + "/" + file + ".h5";

        // handle for file, dataset, and dataspace
        hid_t   fileID, groupID, dsetID, dspID;
        herr_t  h5err; // handle for returned error
        std::vector<hsize_t>     size(1); // for setting data space

        // create a HDF5 file based on the input file name
        fileID = H5Fcreate(extFile.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

        // loop through all fields
        for(unsigned int comp=0; comp<5; ++comp)
        {
            // skip this iteration if no w-velocity
            if ((comp == 2) && (dim == 2)) continue; 

            // the group name of this field
            std::string     gName("/" + fd2str[Field(comp)]);

            // create a group for pressure grid
            groupID = H5Gcreate(
                    fileID, gName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

            // loop through all direction; even in 2D, some post-processing file
            // formats still require z-coordinates, so we aleays output z-coord
            for(unsigned int dir=0; dir<3; ++dir)
            {
                // pre-create a space for the data
                size[0] = n[comp][dir];
                dspID = H5Screate_simple(1, size.data(), nullptr);

                // create a hdf5 data set
                std::string     dsetName(gName + "/" + dir2str[Dir(dir)]);
                dsetID = H5Dcreate2(groupID, dsetName.c_str(), H5T_NATIVE_DOUBLE,
                        dspID, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

                // write the dataset to the h5 file
                h5err = H5Dwrite(dsetID, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
                        H5P_DEFAULT, coord[comp][dir].data());

                // close the data space handle
                h5err = H5Sclose(dspID);
                // close the data set handle
                h5err = H5Dclose(dsetID);
            }

            // close the group for pressure grid
            h5err = H5Gclose(groupID);
        }

        // close the file
        h5err = H5Fclose(fileID); 
    }

    ierr = MPI_Barrier(*comm); CHKERRQ(ierr);

    PetscFunctionReturn(0);
# endif
}


PetscErrorCode CartesianMesh::generateXDMF(
        const std::string &xml, const std::string &file) const
{
    PetscFunctionBeginUser;

    PetscErrorCode  ierr;

    if (mpiRank == 0)
    {
        std::cout << "creating " << xml+".xmf" 
            << " for grid file " << file << " ... ";

        std::ofstream       fs(xml+".xmf");
        fs << "<?xml version=\"1.0\" ?>" << std::endl;
        fs << "<Xdmf Version=\"2.0\">" << std::endl;
        fs << "\t<Domain>" << std::endl;
        fs << "\t\t<Grid Name=\"Vertex Grid\" GridType=\"Uniform\">" << std::endl;

        fs << "\t\t\t<Topology TopologyType=\"3DRectMesh\" NumberOfElements=\"";
        fs << n[4][2] << " " << n[4][1] << " " << n[4][0] << "\"/>" << std::endl;

        fs << "\t\t\t<Geometry GeometryType=\"VxVyVz\">" << std::endl;

        fs << "\t\t\t\t<DataItem Dimensions=\"" << n[4][0] << "\" ";
        fs << "NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
        fs << "\t\t\t\t\t" << file << ":/vertex/x" << std::endl;
        fs << "\t\t\t\t</DataItem>" << std::endl;

        fs << "\t\t\t\t<DataItem Dimensions=\"" << n[4][1] << "\" ";
        fs << "NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
        fs << "\t\t\t\t\t" << file << ":/vertex/y" << std::endl;
        fs << "\t\t\t\t</DataItem>" << std::endl;

        fs << "\t\t\t\t<DataItem Dimensions=\"" << n[4][2] << "\" ";
        fs << "NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
        fs << "\t\t\t\t\t" << file << ":/vertex/z" << std::endl;
        fs << "\t\t\t\t</DataItem>" << std::endl;

        fs << "\t\t\t</Geometry>" << std::endl;
        fs << "\t\t</Grid>" << std::endl;
        fs << "\t</Domain>" << std::endl;
        fs << "</Xdmf>" << std::endl;

        std::cout << "done." << std::endl;
    }

    ierr = MPI_Barrier(*comm); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


PetscErrorCode CartesianMesh::initDMDA()
{
    PetscFunctionBeginUser;

    PetscErrorCode  ierr;

    ierr = createLambdaPack(); CHKERRQ(ierr);
    ierr = createVelocityPack(); CHKERRQ(ierr);

    // total number of local velocity points
    qNLocal = m[0][0] * m[0][1] * m[0][2] + 
        m[1][0] * m[1][1] * m[1][2] + m[2][0] * m[2][1] * m[2][2];

    // total number of local pressure points
    // TODO: remember to change this when we have bodies in the future
    lambdaNLocal = m[3][0] * m[3][1] * m[3][2];

    ierr = createMapping(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


PetscErrorCode CartesianMesh::createSingleDMDA(const PetscInt &i)
{
    using namespace std;

    PetscFunctionBeginUser;

    PetscErrorCode  ierr;

    map<BCType, DMBoundaryType>  petscBC {{PERIODIC, DM_BOUNDARY_PERIODIC},
        {NEUMANN, DM_BOUNDARY_GHOSTED}, {DIRICHLET, DM_BOUNDARY_GHOSTED},
        {CONVECTIVE, DM_BOUNDARY_GHOSTED}};

    DMDALocalInfo   lclInfo;

    switch (dim)
    {
        case 2:
            ierr = DMDACreate2d(*comm,
                    petscBC[(*bcInfo)[XMINUS][u].type],
                    petscBC[(*bcInfo)[YMINUS][v].type], 
                    DMDA_STENCIL_STAR, 
                    n[i][0], n[i][1], nProc[0], nProc[1],
                    1, 1, nullptr, nullptr, &da[i]); CHKERRQ(ierr);
            break;
        case 3:
            ierr = DMDACreate3d(*comm, 
                    petscBC[(*bcInfo)[XMINUS][u].type],
                    petscBC[(*bcInfo)[YMINUS][v].type], 
                    petscBC[(*bcInfo)[ZMINUS][w].type], 
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


PetscErrorCode CartesianMesh::createLambdaPack()
{
    using namespace std;

    PetscFunctionBeginUser;

    PetscErrorCode  ierr;

    ierr = createSingleDMDA(3); CHKERRQ(ierr);

    ierr = DMDAGetInfo(da[3], nullptr, nullptr, nullptr, nullptr,
           &nProc[0], &nProc[1], &nProc[2], 
           nullptr, nullptr, nullptr, nullptr, nullptr, nullptr); CHKERRQ(ierr);

    ierr = DMCompositeCreate(*comm, &lambdaPack); CHKERRQ(ierr);
    ierr = DMCompositeAddDM(lambdaPack, da[3]); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


PetscErrorCode CartesianMesh::createVelocityPack()
{
    using namespace std;

    PetscFunctionBeginUser;

    PetscErrorCode  ierr;

    ierr = DMCompositeCreate(*comm, &qPack); CHKERRQ(ierr);

    for(unsigned int i=0; i<dim; ++i)
    {
        ierr = createSingleDMDA(i); CHKERRQ(ierr);
        ierr = DMCompositeAddDM(qPack, da[i]); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}


PetscErrorCode CartesianMesh::createMapping()
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    ISLocalToGlobalMapping      *mapping;

    ierr = DMCompositeGetISLocalToGlobalMappings(qPack, &mapping); CHKERRQ(ierr);

    // not a hard copy, so do not free the memory space 
    // before we destroy this CartesianMesh instance ...
    qMapping.assign(mapping, mapping+dim);


    // there is a bug (or they do that with a purpose I don't know) in the 
    // mapping returned by DMCompositeLocalToGlobalMappings. That is, the 
    // global indices of ghost points are not -1. Their values are the offset
    // instead. So we have to change these values to -1 manually.
    for(PetscInt c=1; c<dim; ++c)
    {
        const PetscInt      *cArry; // PETSc only return const array ...
        PetscInt            *arry, n, smallest;

        ierr = ISLocalToGlobalMappingGetSize(qMapping[c], &n); CHKERRQ(ierr);
        ierr = ISLocalToGlobalMappingGetIndices(qMapping[c], &cArry); CHKERRQ(ierr);

        // cast to normal ptr that can be modified
        arry = const_cast<PetscInt*>(cArry);

        // Currently, PETSc seems to use the same ghost index on all
        // processes, so we don't have to communicate with other processes
        // to find out the real minimum across them. But we should
        // keep an eye on it to avoid PETSc changes that in the future.
        smallest = *std::min_element(arry, arry+n);

        // if an entry has the min value, it should be a ghost point
        for(PetscInt i=0; i<n; ++i) if (arry[i] == smallest) arry[i] = -1;

        ierr = ISLocalToGlobalMappingRestoreIndices(qMapping[c], &cArry); CHKERRQ(ierr);
    }

    ierr = DMCompositeGetISLocalToGlobalMappings(lambdaPack, &mapping); CHKERRQ(ierr);

    // not a hard copy
    // TODO: now we don't consider bodies. In the future when we have bodies,
    // we may want to modify this.
    lambdaMapping.assign(mapping, mapping+1);

    PetscFunctionReturn(0);
}


std::ostream &operator<< (std::ostream &os, const CartesianMesh &mesh)
{
    os << mesh.info;
    return os;
}
