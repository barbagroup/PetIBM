/***************************************************************************//**
 * \file CartesianMesh.h
 * \author Anush Krishnan (anus@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Definition of the class `CartesianMesh`.
 */


# pragma once

// here goes C++ STL
# include <string>
# include <vector>
# include <map>
# include <memory>
# include <functional>

// here goes PETSc headers
# include <petscsys.h>
# include <petscao.h>
# include <petscdmda.h>
# include <petscdmcomposite.h>

// here goes YAML header
# include <yaml-cpp/yaml.h>

// here goes headers from our PetIBM
# include "types.h"


namespace petibm
{
namespace utilities
{

// TODO: if we can get the column index of ANY velocity point, should we still
//       use ISLocalToGlobalMapping for creating operators?
/** \brief class of composite Cartesian meshes of different field. */
class CartesianMesh
{
public:

    /** \brief dimension. */
    PetscInt                dim = -1;

    /** \brief minimum coordinates of boundaries in all directions. */
    types::RealVec1D        min;

    /** \brief maximum coordinates of boundaries in all directions. */
    types::RealVec1D        max;

    /** \brief total number of points of all fields and in all directions. */
    types::IntVec2D         n;

    /** \brief coordinates of mesh points of all fields and in all directions. */
    types::GhostedVec3D     coord;

    /** \brief point spacing of mesh points of all fields and in all directions. */
    types::GhostedVec3D     dL;

    /** \brief total number of velocity points. */
    PetscInt                UN;

    /** \brief total number of pressure points. */
    PetscInt                pN;

    /** \brief a string for printing information. */
    std::string             info;



    // PETSc stuffs
    /** \brief a vector of DMs of all fields. */
    std::vector<DM>         da;

    /** \brief number of processes in all directions. */
    types::IntVec1D         nProc;

    /** \brief the beginging index of all fields in all directions of this process. */
    types::IntVec2D         bg;

    /** \brief the ending index of all fields in all directions of this process. */
    types::IntVec2D         ed;

    /** \brief the number of points of all fields in all directions of this process. */
    types::IntVec2D         m;

    /** \brief total number of velocity points local to this process. */
    PetscInt                UNLocal;

    /** \brief total number of pressure points local to this process. */
    PetscInt                pNLocal;

    /** \brief DMComposte of velocity DMs. */
    DM                      UPack;

    /** \brief mapping between local unpacked and global packed indices of veloicity.*/
    std::vector<ISLocalToGlobalMapping>     UMapping;

    /** \brief mapping between local and global indices of pressure.*/
    ISLocalToGlobalMapping                  pMapping;



    // MPI stuffs
    /** \brief reference to communicator. */
    std::shared_ptr<const MPI_Comm>         comm;

    /** \brief total number of processes. */
    PetscMPIInt                             mpiSize;

    /** \brief rank of this process. */
    PetscMPIInt                             mpiRank;

    /** \brief reference to BC information. */
    std::shared_ptr<types::BCInfoHolder>    bcInfo;



    /** \brief default constructor. */
    CartesianMesh();


    /**
     * \brief constructor.
     *
     * \param world MPI communicator.
     * \param node a YAML node containing setting of Cartesian mesh.
     * \param bcInfo a BCInfoHolder instance provided by FlowDescription instance.
     * \param type output type.
     */
    CartesianMesh(const MPI_Comm &world, 
            const YAML::Node &node, types::BCInfoHolder &bcInfo, 
            const types::OutputType &type=types::HDF5);


    /** \brief default destructor. */
    ~CartesianMesh();



    /**
     * \brief initialization.
     *
     * \param world MPI communicator.
     * \param node a YAML node containing setting of Cartesian mesh.
     * \param bcInfo a BCInfoHolder instance provided by FlowDescription instance.
     * \param type output type.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode init(const MPI_Comm &world,
            const YAML::Node &meshNode, types::BCInfoHolder &bcInfo, 
            const types::OutputType &type=types::HDF5);


    /**
     * \brief print information.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode printInfo() const;


    /**
     * \brief get the local index of a velocity point by providing MatStencil.
     *
     * Note that the stencil must belong to this process.
     *
     * \param f target field (u=0, v=1, w=2).
     * \param s MatStencil of target velocity point.
     * \param idx returned local index.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode getLocalIndex(
            const PetscInt &f, const MatStencil &s, PetscInt &idx) const;


    /**
     * \brief get the local index of a velocity point by i, j, k.
     *
     * Note that the stencil must belong to this process. For 2D mesh, the
     * value of k-index will be ignored.
     *
     * \param f target field (u=0, v=1, w=2).
     * \param i i-index.
     * \param j j-index.
     * \param k k-index.
     * \param idx returned local index.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode getLocalIndex(const PetscInt &f, 
            const PetscInt &i, const PetscInt &j, const PetscInt &k, 
            PetscInt &idx) const;


    /**
     * \brief get the natural index of a velocity point by providing MatStencil.
     *
     * \param f target field (u=0, v=1, w=2).
     * \param s MatStencil of target velocity point.
     * \param idx returned natural index.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode getNaturalIndex(
            const PetscInt &f, const MatStencil &s, PetscInt &idx) const;


    /**
     * \brief get the natural index of a velocity point by i, j, k.
     *
     * For 2D mesh, the value of k-index will be ignored.
     *
     * \param f target field (u=0, v=1, w=2).
     * \param i i-index.
     * \param j j-index.
     * \param k k-index.
     * \param idx returned natural index.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode getNaturalIndex(const PetscInt &f, 
            const PetscInt &i, const PetscInt &j, const PetscInt &k, 
            PetscInt &idx) const;


    /**
     * \brief get the global index of a velocity point in un-packed DM by 
     *        providing MatStencil.
     *
     * \param f target field (u=0, v=1, w=2).
     * \param s MatStencil of target velocity point.
     * \param idx returned global index in un-packed DM.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode getGlobalIndex(
            const PetscInt &f, const MatStencil &s, PetscInt &idx) const;


    /**
     * \brief get the global index of a velocity point in un-packed DM by 
     *        i, j, k.
     *
     * For 2D mesh, the value of k-index will be ignored.
     *
     * \param f target field (u=0, v=1, w=2).
     * \param i i-index.
     * \param j j-index.
     * \param k k-index.
     * \param idx returned global index in un-packed DM.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode getGlobalIndex(const PetscInt &f, 
            const PetscInt &i, const PetscInt &j, const PetscInt &k, 
            PetscInt &idx) const;


    /**
     * \brief get the global index of a velocity point in packed DM by 
     *        providing MatStencil.
     *
     * \param f target field (u=0, v=1, w=2).
     * \param s MatStencil of target velocity point.
     * \param idx returned global index in packed DM.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode getPackedGlobalIndex(
            const PetscInt &f, const MatStencil &s, PetscInt &idx) const;


    /**
     * \brief get the global index of a velocity point in packed DM by 
     *        i, j, k.
     *
     * For 2D mesh, the value of k-index will be ignored.
     *
     * \param f target field (u=0, v=1, w=2).
     * \param i i-index.
     * \param j j-index.
     * \param k k-index.
     * \param idx returned global index in packed DM.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode getPackedGlobalIndex(const PetscInt &f, 
            const PetscInt &i, const PetscInt &j, const PetscInt &k, 
            PetscInt &idx) const;


    /**
     * \brief set output format.
     *
     * Through this function, programs can change the type of output after 
     * initialization.
     *
     * \param type output type.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode setOutputFormat(const types::OutputType &type);


    /**
     * \brief create XDMF file for grid only.
     *
     * \param xml the file name of XDMF file, path included.
     * \param file the file name of HDF5 file, path included.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode generateXDMF(
            const std::string &xml, const std::string &file) const;


    /** \brief a function for writing output data.
     *
     * This is linked to a real writing function based on the type of output. 
     * The first parameter is the path of folder for the output file. The second
     * parameter is the file name without extension. The extension will be 
     * determined by the type of output.
     */
    std::function<PetscErrorCode(const std::string &, const std::string &)> write;

protected:

    /** \brief the underlying data for mesh point spacing. */
    types::RealVec3D        dLTrue;

    /** \brief the underlying data for mesh point coordinates. */
    types::RealVec3D        coordTrue;

    /** \brief references to underlying AO objects of velocity DMs. */
    std::vector<AO>         ao;

    /** \brief number of local velocity points (without ghosts) for all 
     *         processes and all velocity fields. */
    types::IntVec2D         UNLocalAllProcs;

    /** \brief number of local packed velocity points (without ghost) for 
     *         all processes. */
    types::IntVec1D         UPackNLocalAllProcs;

    /** \brief offsets of velocity points in un-packed DMs for each field. */
    types::IntVec2D         offsetsAllProcs;

    /** \brief offsets of packed velocity points in packed DM. */
    types::IntVec1D         offsetsPackAllProcs;

    /**
     * \brief create vertex information.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode createVertexMesh();


    /**
     * \brief create pressure mesh information.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode createPressureMesh();


    /**
     * \brief create velocity mesh information.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode createVelocityMesh();


    /**
     * \brief create a string of information.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode createInfoString();


    /**
     * \brief gather information of parallel distribution and add to info string.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode addLocalInfoString(std::stringstream &ss);


    /**
     * \brief underlying function to write binary mesh data.
     *
     * \param dir the path of the folder for output file.
     * \param file the file name without extension.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode writeBinary(const std::string &dir, const std::string &file);


    /**
     * \brief underlying function to write VTK mesh file.
     *
     * \param dir the path of the folder for output file.
     * \param file the file name without extension.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode writeVTK(const std::string &dir, const std::string &file);


    /**
     * \brief underlying function to write HDF5 mesh data.
     *
     * \param dir the path of the folder for output file.
     * \param file the file name without extension.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode writeHDF5(const std::string &dir, const std::string &file);



    /**
     * \brief initialize DMDAs.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode initDMDA();


    /**
     * \brief function for creating a single DMDA.
     *
     * \param i the index of the targeting field (0 ~ 4 represents u, v, w, 
     *          pressure, and vertex respectively.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode createSingleDMDA(const PetscInt &i);


    /**
     * \brief create DMDA for pressure.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode createPressureDMDA();


    /**
     * \brief create DMDAs for velovity fields and make a DMComposite.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode createVelocityPack();


    /**
     * \brief create ISLocalToGlobalMapping for fileds.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode createMapping();

}; // CartesianMesh


/**
 * \brief an I/O function for using standard output stream.
 *
 * \param os output stream.
 * \param mesh an instance of CartesainMesh.
 *
 * \return output stream.
 */
std::ostream &operator<< (std::ostream &os, const CartesianMesh &mesh);

} // end of namespace utilities
} // end of namespace petibm
