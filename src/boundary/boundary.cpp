/*
 * boundary.cpp
 * Copyright (C) 2017 Pi-Yueh Chuang <pychuang@gwu.edu>
 *
 * Distributed under terms of the MIT license.
 */

# include <petibm/boundary.h>
# include <petibm/boundarysimple.h>


namespace petibm
{
namespace boundary
{

PetscErrorCode createBoundary(
        const type::Mesh &mesh, const YAML::Node &node,
        type::Boundary &boundary)
{
    PetscFunctionBeginUser;
    
    boundary = std::make_shared<BoundarySimple>(mesh, node);
    
    PetscFunctionReturn(0);
}

}
}
