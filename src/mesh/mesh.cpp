/*
 * mesh.cpp
 * Copyright (C) 2017 Pi-Yueh Chuang <pychuang@gwu.edu>
 *
 * Distributed under terms of the MIT license.
 */

# include <petibm/mesh.h>
# include <petibm/cartesianmesh.h>
# include <petibm/io.h>


namespace petibm
{
namespace mesh
{
    PetscErrorCode MeshBase::printInfo() const
    {
        PetscFunctionBeginUser;
        
        PetscErrorCode ierr;
        ierr = io::print(info); CHKERRQ(ierr);
        
        PetscFunctionReturn(0);
    }
    
    
    // factory function to create Mesh.
    PetscErrorCode createMesh(
            const MPI_Comm &comm, const YAML::Node &node, type::Mesh &mesh)
    {
        PetscFunctionBeginUser;
        
        mesh = std::make_shared<CartesianMesh>(comm, node);
        
        PetscFunctionReturn(0);
    }
    
} // end of mesh
} // end of petibm
