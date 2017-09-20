/*
 * mesh.cpp
 * Copyright (C) 2017 Pi-Yueh Chuang <pychuang@gwu.edu>
 *
 * Distributed under terms of the MIT license.
 */

#include <petibm/mesh.h>
#include <petibm/cartesianmesh.h>


namespace petibm
{
namespace mesh
{
    
    // factory function to create Mesh.
    PetscErrorCode createMesh(
            const MPI_Comm &comm, const YAML::Node &node, type::Mesh &mesh)
    {
        PetscFunctionBeginUser;
        
        mesh = type::Mesh(new CartesianMesh(comm, node));
        
        PetscFunctionReturn(0);
    }
    
} // end of mesh
} // end of petibm
