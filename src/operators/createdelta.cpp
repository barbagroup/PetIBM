/**
 * \file createdelta.cpp
 * \brief Definition of functions creating Delta operator.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 */

#include <cmath>

#include <petscmat.h>

#include <petibm/bodypack.h>
#include <petibm/boundary.h>
#include <petibm/delta.h>
#include <petibm/mesh.h>
#include <petibm/singlebody.h>
#include <petibm/type.h>

namespace petibm
{
namespace operators
{
// TODO: it's anti-readable that we mix the use of local and global Lagrangian
// index
// TODO: no pre-allocation for D matrix, this may be inefficient, though it
// works

PetscErrorCode getEulerianNeighbors(
    const type::Mesh &mesh, const PetscInt &dof, const type::IntVec1D &IJK,
    const std::vector<bool> &periodic, const PetscInt &window,
    type::IntVec2D &ijk, type::RealVec2D &xyz);

// implementation of petibm::operators::createDelta
PetscErrorCode createDelta(const type::Mesh &mesh, const type::Boundary &bc,
                           const type::BodyPack &bodies,
                           const delta::DeltaKernel &kernel,
                           const PetscInt &kernelSize,
                           Mat &Op)
{
    PetscFunctionBeginUser;

    PetscErrorCode ierr;

    // get periodic flags and domain sizes
    std::vector<bool> periodic(mesh->dim);  // flags to check periodicity
    for (PetscInt d = 0; d < mesh->dim; ++d)
    {
        // the the x-component of the velocity is periodic in a direction,
        // so are the other components of the velocity
        periodic[d] = (bc->bds[0][d]->type == type::BCType::PERIODIC);
    }

    ierr = MatCreate(mesh->comm, &Op); CHKERRQ(ierr);
    ierr = MatSetSizes(Op, bodies->nLclPts * bodies->dim, mesh->UNLocal,
                       bodies->nPts * bodies->dim, mesh->UN); CHKERRQ(ierr);
    ierr = MatSetFromOptions(Op); CHKERRQ(ierr);
    ierr = MatSetUp(Op); CHKERRQ(ierr);
    ierr = MatSetOption(
        Op, MAT_KEEP_NONZERO_PATTERN, PETSC_FALSE); CHKERRQ(ierr);
    ierr = MatSetOption(
        Op, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE); CHKERRQ(ierr);

    // loop through all bodies
    for (PetscInt bIdx = 0; bIdx < bodies->nBodies; ++bIdx)
    {
        // get an alias of the current body (for code readability)
        const type::SingleBody &body = bodies->bodies[bIdx];

        // get the cell widths of the background mesh
        // (the regularized delta functions work for uniform meshes)
        std::vector<PetscReal> widths(mesh->dim);
        for (PetscInt d = 0; d < mesh->dim; ++d)
        {
            // get the directional width of the Eulerian cell
            // closest to the first Lagrangian point
            widths[d] = mesh->dL[0][d][body->meshIdx[0][d]];
        }

        // loop through all local points of this body
        for (PetscInt iLcl = 0, iGlb = body->bgPt; iLcl < body->nLclPts;
             iLcl++, iGlb++)
        {
            // get alias of coordinates and background index of current point
            const type::IntVec1D &IJK = body->meshIdx[iLcl];
            const type::RealVec1D &XYZ = body->coords[iGlb];

            // loop through all directions
            for (PetscInt dof = 0; dof < body->dim; ++dof)
            {
                type::IntVec1D cols;
                type::RealVec1D vals;

                // get packed row index
                PetscInt row;
                ierr = bodies->getPackedGlobalIndex(
                    bIdx, iGlb, dof, row); CHKERRQ(ierr);

                type::IntVec2D ijk;
                type::RealVec2D coords;
                ierr = getEulerianNeighbors(
                    mesh, dof, IJK, periodic, kernelSize, ijk, coords); CHKERRQ(ierr);

                type::RealVec1D xyz(mesh->dim, 0.0);
                if (mesh->dim == 3)
                {
                    for (unsigned int k = 0; k < ijk[2].size(); ++k)
                    {
                        xyz[2] = coords[2][k];
                        for (unsigned int j = 0; j < ijk[1].size(); ++j)
                        {
                            xyz[1] = coords[1][j];
                            for (unsigned int i = 0; i < ijk[0].size(); ++i)
                            {
                                xyz[0] = coords[0][i];
                                // get packed column index
                                PetscInt col;
                                ierr = mesh->getPackedGlobalIndex(
                                    dof, ijk[0][i], ijk[1][j], ijk[2][k],
                                    col); CHKERRQ(ierr);
                                
                                PetscReal val;
                                val = delta::delta(
                                    XYZ, xyz, widths, kernel); CHKERRQ(ierr);
                                
                                cols.push_back(col);
                                vals.push_back(val);
                            }
                        }
                    }
                }
                else if (mesh->dim == 2)
                {
                    for (unsigned int j = 0; j < ijk[1].size(); ++j)
                    {
                        xyz[1] = coords[1][j];
                        for (unsigned int i = 0; i < ijk[0].size(); ++i)
                        {
                            xyz[0] = coords[0][i];
                            // get packed column index
                            PetscInt col;
                            ierr = mesh->getPackedGlobalIndex(
                                dof, ijk[0][i], ijk[1][j], 0,
                                col); CHKERRQ(ierr);
                            
                            PetscReal val;
                            val = delta::delta(
                                XYZ, xyz, widths, kernel); CHKERRQ(ierr);
                            
                            cols.push_back(col);
                            vals.push_back(val);
                        }
                    }
                }
                else
                    SETERRQ(mesh->comm, PETSC_ERR_ARG_WRONG,
                            "Only 2D and 3D configurations are supported.\n");

                ierr = MatSetValues(Op, 1, &row,
                                    cols.size(), cols.data(), vals.data(),
                                    INSERT_VALUES); CHKERRQ(ierr);
            }
        }
    }

    ierr = MatAssemblyBegin(Op, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Op, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // createDelta

PetscErrorCode getEulerianNeighbors(
    const type::Mesh &mesh, const PetscInt &dof, const type::IntVec1D &IJK,
    const std::vector<bool> &periodic, const PetscInt &window,
    type::IntVec2D &ijk, type::RealVec2D &xyz)
{
    PetscFunctionBeginUser;

    xyz.resize(mesh->dim);
    ijk.resize(mesh->dim);

    for (PetscInt d = 0; d < mesh->dim; ++d)
    {
        for (PetscInt s = IJK[d] - window; s <= IJK[d] + window; ++s)
        {
            if (s >= 0 && s < mesh->n[dof][d])
            {
                ijk[d].push_back(s);
                xyz[d].push_back(mesh->coord[dof][d][s]);
            }
            else if (periodic[d])
            {
                PetscReal L = mesh->max[d] - mesh->min[d];
                if (s < 0)
                {
                    ijk[d].push_back(s + mesh->n[dof][d]);
                    xyz[d].push_back(mesh->coord[dof][d][s] - L);
                }
                else
                {
                    ijk[d].push_back(s - mesh->n[dof][d]);
                    xyz[d].push_back(mesh->coord[dof][d][s] + L);
                }
            }
        }
    }

    PetscFunctionReturn(0);
}  // getEulerianNeighbors

}  // end of namespace operators
}  // end of namespace petibm
