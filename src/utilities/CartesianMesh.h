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

// here goes PETSc headers
# include <petscsys.h>

// here goes YAML header
# include <yaml-cpp/yaml.h>

// here goes headers from our PetIBM
# include "types.h"


class CartesianMesh
{
public:

    PetscInt                dim = -1;

    types::RealVec1D        bg, ed;

    types::IntVec2D         n = types::IntVec2D(5);

    types::RealVec3D        coord = types::RealVec3D(5);

    types::RealVec3D        dL = types::RealVec3D(5);

    std::string             info;



    CartesianMesh();
    CartesianMesh(const std::string &file, types::BCInfoHolder &bcInfo);
    CartesianMesh(const YAML::Node &node, types::BCInfoHolder &bcInfo);

    ~CartesianMesh();

    PetscErrorCode init(
            const YAML::Node &meshNode, types::BCInfoHolder &bcInfo);

    PetscErrorCode printInfo();

    template <types::OutputType out>
    PetscErrorCode write(const std::string &file, const std::string &xml="");

    PetscErrorCode generateXDMF(const std::string &xml, const std::string &file);

protected:

    PetscErrorCode createVertexMesh();
    PetscErrorCode createPressureMesh();
    PetscErrorCode createVelocityMesh(types::BCInfoHolder &bcInfo);
    PetscErrorCode createInfoString();

}; // CartesianMesh


std::ostream &operator<< (std::ostream &os, const CartesianMesh &mesh);
