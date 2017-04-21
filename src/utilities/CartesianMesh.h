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
# include <petscdmda.h>
# include <petscdmcomposite.h>

// here goes YAML header
# include <yaml-cpp/yaml.h>

// here goes headers from our PetIBM
# include "types.h"


class CartesianMesh
{
public:

    PetscInt                dim = -1;

    types::RealVec1D        min,
                            max;

    types::IntVec2D         n;

    types::RealVec3D        coord;

    types::DeltaLVec        dL;

    std::string             info;



    // PETSc stuffs
    std::vector<DM>         da;

    types::IntVec1D         nProc;

    types::IntVec2D         bg;

    types::IntVec2D         ed;

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
    std::shared_ptr<const MPI_Comm>         comm;

    /** \brief total number of processes. */
    PetscMPIInt                             mpiSize;

    /** \brief rank of this process. */
    PetscMPIInt                             mpiRank;

    // a reference to BC information
    std::shared_ptr<types::BCInfoHolder>    bcInfo;


    // cunstructors
    CartesianMesh();

    CartesianMesh(const MPI_Comm &world, 
            const YAML::Node &node, types::BCInfoHolder &bcInfo, 
            const types::OutputType &type=types::HDF5);

    ~CartesianMesh();

    // real initialization function
    PetscErrorCode init(const MPI_Comm &world,
            const YAML::Node &meshNode, types::BCInfoHolder &bcInfo, 
            const types::OutputType &type=types::HDF5);

    // here goes the part regarding PETSc DMDA objects
    PetscErrorCode initDMDA();

    // print info with PetscPrintf
    PetscErrorCode printInfo() const;

    // to change the OutputType in after the instance was initialized
    PetscErrorCode setOutputFormat(const types::OutputType &type);

    // a function to create XDMF file for HDF5 grid file
    PetscErrorCode generateXDMF(
            const std::string &xml, const std::string &file) const;

    // a function for writing out the grid based on the OutputType specified
    std::function<PetscErrorCode(const std::string &, const std::string &)> write;

protected:

    types::RealVec3D        dLTrue;

    PetscErrorCode createVertexMesh();
    PetscErrorCode createPressureMesh();
    PetscErrorCode createVelocityMesh();
    PetscErrorCode createInfoString();
    PetscErrorCode addLocalInfoString(std::stringstream &ss);

    PetscErrorCode writeBinary(const std::string &dir, const std::string &file);
    PetscErrorCode writeVTK(const std::string &dir, const std::string &file);
    PetscErrorCode writeHDF5(const std::string &dir, const std::string &file);

    PetscErrorCode createSingleDMDA(const PetscInt &i);
    PetscErrorCode createPressureDMDA();
    PetscErrorCode createVelocityPack();
    PetscErrorCode createMapping();

}; // CartesianMesh


std::ostream &operator<< (std::ostream &os, const CartesianMesh &mesh);
