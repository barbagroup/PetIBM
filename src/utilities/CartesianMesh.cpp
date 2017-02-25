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
CartesianMesh::CartesianMesh() {};


CartesianMesh::CartesianMesh(const std::string &file, BCInfoHolder &bc, const OutputType &type)
{
    YAML::Node  node = YAML::LoadFile(file);
    
    CartesianMesh(node, bc);
}


CartesianMesh::CartesianMesh(const YAML::Node &node, BCInfoHolder &bc, const OutputType &type)
{
    if (node["cartesianMesh"].IsDefined())
        init(node["cartesianMesh"], bc, type);
    else
        init(node, bc, type);
}


CartesianMesh::~CartesianMesh() {};


PetscErrorCode CartesianMesh::init(
        const YAML::Node &meshNode, BCInfoHolder &bc, const OutputType &type)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    ierr = PetscPrintf(PETSC_COMM_WORLD, "Creating a CartesianMesh object ... ");
    CHKERRQ(ierr);

    // index 3 represent pressure mesh; min & max always represent pressure mesh
    ierr = parser::parseMesh(meshNode, dim, min, max, n[3], dL[3]); CHKERRQ(ierr);

    // store the address of bcInfo
    // note: this is a bad practice; shared_ptr is not for stack variables!!
    bcInfo = std::shared_ptr<BCInfoHolder>(&bc, [](BCInfoHolder*){});

    ierr = createPressureMesh(); CHKERRQ(ierr);
    ierr = createVertexMesh(); CHKERRQ(ierr);
    ierr = createVelocityMesh(); CHKERRQ(ierr);
    ierr = createInfoString(); CHKERRQ(ierr);

    ierr = setOutputFormat(type); CHKERRQ(ierr);

    ierr = initDMDA(); CHKERRQ(ierr);

    ierr = PetscPrintf(PETSC_COMM_WORLD, "done.\n"); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


PetscErrorCode CartesianMesh::setOutputFormat(const OutputType &type)
{
    PetscFunctionBeginUser;

    using namespace std::placeholders;

    switch (type)
    {
        case OutputType::Binary:
            write = std::bind(&CartesianMesh::writeBinary, this, _1);
            break;
        case OutputType::VTK:
            write = std::bind(&CartesianMesh::writeVTK, this, _1);
            break;
        case OutputType::HDF5:
            write = std::bind(&CartesianMesh::writeHDF5, this, _1);
            break;
    }

    PetscFunctionReturn(0);
}


PetscErrorCode CartesianMesh::createPressureMesh()
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    // note: no matter it's 2D or 3D, we always create z coordinates
    
    // initialization; other property should be already created in parseMesh
    // note: index 3 means pressure mesh
    coord[3] = {RealVec1D(1, 0.0), RealVec1D(1, 0.0), RealVec1D(1, 0.0)};

    // loop through all axes to get the coordinates of pressure points
    for(unsigned int i=0; i<dim; ++i)
    {
        coord[3][i].resize(n[3][i]);

        // the following three lines calculate the coord of pressure points
        std::partial_sum(dL[3][i].begin(), dL[3][i].end(), coord[3][i].begin());

        auto f = [i, this] (double &a, double &b) -> double 
                    { return a + this->min[i] - 0.5 * b; };

        std::transform(coord[3][i].begin(), coord[3][i].end(), 
                dL[3][i].begin(), coord[3][i].begin(), f);
    }

    PetscFunctionReturn(0);
}


PetscErrorCode CartesianMesh::createVertexMesh()
{
    PetscFunctionBeginUser;

    // initialization; note: index 4 means vertex mesh
    n[4] = {1, 1, 1};
    coord[4] = {RealVec1D(1, 0.0), RealVec1D(1, 0.0), RealVec1D(1, 0.0)};
    dL[4] = {RealVec1D(1, 0.0), RealVec1D(1, 0.0), RealVec1D(1, 0.0)};

    // loop through all axes to get the coordinates of mesh vertexes
    // note: index 3 means pressure mesh; 4 means vertexes
    for(unsigned int i=0; i<dim; ++i)
    {
        // number of vertexes is one more than pressure cells
        n[4][i] = n[3][i] + 1;
        coord[4][i].resize(n[4][i]);

        // create coordinates of vertexes
        std::partial_sum(dL[3][i].begin(), dL[3][i].end(), coord[4][i].begin()+1);
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
        // initialization
        n[comp] = {1, 1, 1};
        coord[comp] = {RealVec1D(1, 0.0), RealVec1D(1, 0.0), RealVec1D(1, 0.0)};
        dL[comp] = {RealVec1D(1, 0.0), RealVec1D(1, 0.0), RealVec1D(1, 0.0)};

        // loop through each direction, i.e., x, y, z axes
        for(unsigned int dir=0; dir<dim; ++dir)
        {
            // when the direction corresponding to velocity component
            if (dir == comp)
            {
                // there will no velocity point on boundary (it ghosted)
                n[comp][dir] = n[3][dir] - 1;

                // coordinates will match that of the vertex
                coord[comp][dir] = types::RealVec1D(
                        coord[4][dir].begin()+1, coord[4][dir].end()-1);

                // cell size will be a half of the sum of adjacent pressure cell
                dL[comp][dir].resize(n[comp][dir]);

                auto f = [] (const PetscReal &x, const PetscReal &y)
                    -> PetscReal { return 0.5 * (x + y); };

                // note: adjacent can not automatically ignore the first element
                std::adjacent_difference(dL[3][dir].begin()+1, dL[3][dir].end(),
                        dL[comp][dir].begin(), f);

                // so we have to handle the first cell manually
                dL[comp][dir][0] = 0.5 * (dL[3][dir][0] + dL[3][dir][1]);

                // if BC in this direction is periodic, we append one more 
                // velocity point at the end
                if ((*bcInfo)[BCLoc(dir*2)][Field(comp)].type == BCType::PERIODIC)
                {
                    n[comp][dir] += 1;
                    coord[comp][dir].push_back(coord[4][dir][n[comp][dir]]);
                    dL[comp][dir].push_back(
                            f(dL[3][dir][0], dL[3][dir][n[comp][dir]-1]));
                }
            }
            // other directions
            else
            {
                n[comp][dir] = n[3][dir];

                // coordinate and cel size will match that of pressure point
                coord[comp][dir] = coord[3][dir];
                dL[comp][dir] = dL[3][dir];
            }
        }
    }

    PetscFunctionReturn(0);
}


PetscErrorCode CartesianMesh::createInfoString()
{
    PetscFunctionBeginUser;

    std::stringstream       ss;

    ss << "Cartesian Staggered Grid:" << std::endl;

    ss << "\tDimension: " << dim << std::endl;
    ss << "\tDomain range: " 
        << min[0] << " -> " << max[0] << "; "
        << min[1] << " -> " << max[3];
   
    if (dim == 3)
        ss << "; " << min[2] << " -> " << max[2];
    ss << std::endl;

    ss << "\tNumber of pressure cells (Nx x Ny" << ((dim==2)? "" : " x Nz") << "): ";
    ss << n[3][0];
    for(unsigned int dir=1; dir<dim; ++dir) ss << " x " << n[3][dir];
    ss << std::endl;

    for(unsigned int comp=0; comp<dim; ++comp)
    {
        ss << "\tNumber of " << fd2str[Field(comp)] 
            << "-velocity cells (Nx x Ny" << ((dim==2)? "" : " x Nz") << "): ";
        ss << n[comp][0];
        for(unsigned int dir=1; dir<dim; ++dir) ss << " x " << n[comp][dir];
        ss << std::endl;
    }

    info = ss.str();

    PetscFunctionReturn(0);
}


PetscErrorCode CartesianMesh::printInfo()
{
    PetscFunctionBeginUser;
    
    PetscErrorCode ierr;
    ierr = PetscPrintf(PETSC_COMM_WORLD, info.c_str()); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // printInfo


PetscErrorCode CartesianMesh::writeBinary(const std::string &file)
{
    PetscErrorCode ierr;

    PetscMPIInt rank;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    if (rank == 0)
    {
        std::ofstream   streamFile(file+".txt");

        std::cout << "Writing grid to TXT file: " << file+".txt" << " ... ";

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

        std::cout << "done." << std::endl;
    }

    ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // write


PetscErrorCode CartesianMesh::writeVTK(const std::string &file)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    // get the rank
    PetscMPIInt         rank;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    // only the master process will write file
    if (rank == 0)
    {
        std::ofstream   fs(file+".vtk");

        std::cout << "Writing grid to VTK file: " << file+".vtk" << " ... ";

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

        std::cout << "done." << std::endl;
    }

    // all processes must to wait the master process
    ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


PetscErrorCode CartesianMesh::writeHDF5(const std::string &file)
{
    PetscFunctionBeginUser;

# ifndef PETSC_HAVE_HDF5
    SETERRQ(PETSC_COMM_WORLD, 56, "Seems the PETSc was not compiled with HDF5.");
# else

    PetscErrorCode  ierr;
    PetscMPIInt     rank;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    if (rank == 0)
    {
        std::string     extFile = file + ".h5";
        std::cout << "Writing grid to HDF5 file: " << extFile+".h5" << " ... ";

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

        std::cout << "done." << std::endl;
    }

    ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);

    PetscFunctionReturn(0);
# endif
}


PetscErrorCode CartesianMesh::generateXDMF(
        const std::string &xml, const std::string &file)
{
    PetscFunctionBeginUser;

    PetscErrorCode  ierr;

    PetscMPIInt     rank;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    if (rank == 0)
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

    ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


PetscErrorCode CartesianMesh::initDMDA()
{
    PetscFunctionBeginUser;

    PetscErrorCode  ierr;

    ierr = createLambdaPack(); CHKERRQ(ierr);
    ierr = createVelocityPack(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


PetscErrorCode CartesianMesh::createSingleDMDA(const PetscInt &i)
{
    using namespace types;
    using namespace std;

    PetscFunctionBeginUser;

    PetscErrorCode  ierr;

    map<BCType, DMBoundaryType>  petscBC {{types::PERIODIC, DM_BOUNDARY_PERIODIC},
        {types::NEUMANN, DM_BOUNDARY_GHOSTED}, {types::DIRICHLET, DM_BOUNDARY_GHOSTED}};

    DMDALocalInfo   lclInfo;

    switch (dim)
    {
        case 2:
            ierr = DMDACreate2d(PETSC_COMM_WORLD,
                    petscBC[(*bcInfo)[XMINUS][u].type],
                    petscBC[(*bcInfo)[XMINUS][v].type], 
                    DMDA_STENCIL_STAR, 
                    n[i][0], n[i][1], PETSC_DECIDE, PETSC_DECIDE,
                    1, 1, nullptr, nullptr, &da[i]); CHKERRQ(ierr);
            break;
        case 3:
            ierr = DMDACreate3d(PETSC_COMM_WORLD, 
                    petscBC[(*bcInfo)[XMINUS][u].type],
                    petscBC[(*bcInfo)[YMINUS][v].type], 
                    petscBC[(*bcInfo)[ZMINUS][w].type], 
                    DMDA_STENCIL_STAR, n[i][0], n[i][1], n[i][2], 
                    PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, 
                    1, 1, nullptr, nullptr, nullptr, &da[i]); CHKERRQ(ierr);
            break;
    }

    ierr = DMDAGetLocalInfo(da[i], &lclInfo); CHKERRQ(ierr);

    bg[i].resize(3); ed[i].resize(3); m[i].resize(3);

    bg[i][0] = lclInfo.xs; bg[i][1] = lclInfo.ys; bg[i][2] = lclInfo.zs;

    m[i][0] = lclInfo.xm; m[i][1] = lclInfo.ym; m[i][2] = lclInfo.zm;

    for(unsigned int comp=0; comp<3; ++comp) 
        ed[i][comp] = bg[i][comp] + m[i][comp];

    PetscFunctionReturn(0);
}


PetscErrorCode CartesianMesh::createLambdaPack()
{
    using namespace types;
    using namespace std;

    PetscFunctionBeginUser;

    PetscErrorCode  ierr;

    ierr = createSingleDMDA(3); CHKERRQ(ierr);

    ierr = DMCompositeCreate(PETSC_COMM_WORLD, &lambdaPack); CHKERRQ(ierr);
    ierr = DMCompositeAddDM(lambdaPack, da[3]); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


PetscErrorCode CartesianMesh::createVelocityPack()
{
    using namespace types;
    using namespace std;

    PetscFunctionBeginUser;

    PetscErrorCode  ierr;

    ierr = DMCompositeCreate(PETSC_COMM_WORLD, &qPack); CHKERRQ(ierr);

    for(unsigned int i=0; i<dim; ++i)
    {
        ierr = createSingleDMDA(i); CHKERRQ(ierr);
        ierr = DMCompositeAddDM(qPack, da[i]); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}


std::ostream &operator<< (std::ostream &os, const CartesianMesh &mesh)
{
    os << mesh.info;
    return os;
}
