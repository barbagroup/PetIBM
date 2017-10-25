/***************************************************************************//**
 * \file createdelta.cpp
 * \author Anush Krishnan (anus@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Definition of functions regarding to creating Delta operators.
 */


// here goes PETSc headers
# include <petscmat.h>

// here goes headers from our PetIBM
# include <petibm/mesh.h>
# include <petibm/boundary.h>
# include <petibm/bodypack.h>
# include <petibm/singlebody.h>
# include <petibm/type.h>
# include <petibm/delta.h>


namespace petibm
{
namespace operators
{

// TODO: it's anti-readiable that we mix the use of local and global Lagrangian index
// TODO: no pre-allocation for D matrix, this may be inefficient, though it works


PetscErrorCode getWindowAndDistance(const PetscInt &dim,
        const type::IntVec1D &n,
        const type::GhostedVec2D &coords,
        const std::vector<bool> &periodic,
        const type::RealVec1D &L,
        const type::IntVec1D &IJK,
        const type::RealVec1D &XYZ,
        type::IntVec2D &targets,
        type::RealVec2D &targetdLs);


/** \copydoc createDelta(const CartesianMesh &, const BodyPack &, Mat &). */
PetscErrorCode createDelta(const type::Mesh &mesh, const type::Boundary &bc,
                           const type::BodyPack &bodies, Mat &D)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    // flags to see if BC is periodic
    std::vector<bool>   periodic(mesh->dim);

    // domain sizes for periodic cases
    type::RealVec1D     L(mesh->dim);

    // get periodic flags and domain sizes
    for(PetscInt d=0; d<mesh->dim; ++d)
    {
        // if u-velocity at a direction is periodic, so are other velocity fields
        periodic[d] = (bc->bds[0][d]->type == type::BCType::PERIODIC);
        
        L[d] = mesh->max[d] - mesh->min[d];
    }


    ierr = MatCreate(mesh->comm, &D); CHKERRQ(ierr);
    ierr = MatSetSizes(D, bodies->nLclPts*bodies->dim, mesh->UNLocal,
            bodies->nPts*bodies->dim, mesh->UN); CHKERRQ(ierr);
    ierr = MatSetFromOptions(D); CHKERRQ(ierr);
    ierr = MatSetUp(D); CHKERRQ(ierr);
    ierr = MatSetOption(D, MAT_KEEP_NONZERO_PATTERN, PETSC_FALSE); CHKERRQ(ierr);
    ierr = MatSetOption(D, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE); CHKERRQ(ierr);


    // loop through all bodies
    for(PetscInt bIdx=0; bIdx<bodies->nBodies; ++bIdx)
    {
        // get an alias of current body for code simplicity
        const type::SingleBody    bd = bodies->bodies[bIdx];

        // loop through all local points of this body
        for(PetscInt iLcl=0, iGlb=bd->bgPt; iLcl<bd->nLclPts; iLcl++, iGlb++)
        {
            // get alias of coordinates and background index of current point
            const type::IntVec1D   &IJK = bd->meshIdx[iLcl];
            const type::RealVec1D  &XYZ = bd->coords[iGlb];

            // loop through all degree of freedom
            for(PetscInt dof=0; dof<bd->dim; ++dof)
            {
                PetscInt            row; // row in packed matrix
                type::IntVec2D      targets; // index of valid velocity points
                type::RealVec2D     targetdLs; // distance to velocity points
                type::IntVec1D      cols;
                type::RealVec1D     values;

                // get packed row index
                ierr = bodies->getPackedGlobalIndex(
                        bIdx, iGlb, dof, row); CHKERRQ(ierr);


                // get indices and distances of valid velocity points
                // valid means the ones that may interact with Lagrangian point
                ierr = getWindowAndDistance(mesh->dim, mesh->n[dof], 
                        mesh->coord[dof], periodic, L, IJK, XYZ, 
                        targets, targetdLs); CHKERRQ(ierr);


                if (mesh->dim == 3)
                {
                    for(unsigned int k=0; k<targets[2].size(); ++k)
                    {
                        const PetscReal &hz = mesh->dL[dof][2][targets[2][k]];

                        for(unsigned int j=0; j<targets[1].size(); ++j)
                        {
                            const PetscReal &hy = mesh->dL[dof][1][targets[1][j]];

                            for(unsigned int i=0; i<targets[0].size(); ++i)
                            {
                                const PetscReal &hx = mesh->dL[dof][0][targets[0][i]];

                                PetscInt    col;
                                PetscReal   value;

                                ierr = mesh->getPackedGlobalIndex(
                                        dof, targets[0][i], targets[1][j], 
                                        targets[0][k], col); CHKERRQ(ierr);

                                value = delta::Roma_et_al( targetdLs[0][i], hx,
                                    targetdLs[1][j], hy, targetdLs[2][k], hz);

                                cols.push_back(col);
                                values.push_back(value);
                            }
                        }
                    }
                }
                else
                {
                    for(unsigned int j=0; j<targets[1].size(); ++j)
                    {
                        const PetscReal &hy = mesh->dL[dof][1][targets[1][j]];

                        for(unsigned int i=0; i<targets[0].size(); ++i)
                        {
                            const PetscReal &hx = mesh->dL[dof][0][targets[0][i]];

                            PetscInt    col;
                            PetscReal   value;

                            ierr = mesh->getPackedGlobalIndex(
                                    dof, targets[0][i], targets[1][j], 
                                    0, col); CHKERRQ(ierr);

                            value = delta::Roma_et_al(
                                    targetdLs[0][i], hx, targetdLs[1][j], hy);

                            cols.push_back(col);
                            values.push_back(value);
                        }
                    }
                }

                ierr = MatSetValues(D, 1, &row, cols.size(), cols.data(), 
                        values.data(), INSERT_VALUES); CHKERRQ(ierr);

            }
        }
    }

    ierr = MatAssemblyBegin(D, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(D, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


PetscErrorCode getWindowAndDistance(const PetscInt &dim,
        const type::IntVec1D &n,
        const type::GhostedVec2D &coords,
        const std::vector<bool> &periodic,
        const type::RealVec1D &L,
        const type::IntVec1D &IJK,
        const type::RealVec1D &XYZ,
        type::IntVec2D &targets,
        type::RealVec2D &targetdLs)
{
    PetscFunctionBeginUser;

    targets.resize(dim);
    targetdLs.resize(dim);

    // loop through all direction
    for(PetscInt dir=0; dir<dim; ++dir)
    {
        // loop through all possible index in this direction (+-2)
        for(PetscInt ijk=IJK[dir]-2; ijk<=IJK[dir]+2; ++ijk)
        {
            if ((ijk >= 0) && (ijk < n[dir]))
            {
                targets[dir].push_back(ijk);
                targetdLs[dir].push_back(XYZ[dir] - coords[dir][ijk]);
            }
            else if (periodic[dir]) // only when periodic, we do something
            {
                if (ijk < 0)
                {
                    targets[dir].push_back(n[dir] + ijk);
                    targetdLs[dir].push_back(XYZ[dir] + L[dir] - 
                            coords[dir][targets[dir].back()]);
                }
                else // imply ijk > # of velocity points
                {
                    targets[dir].push_back(ijk - n[dir]);
                    targetdLs[dir].push_back(XYZ[dir] - L[dir] -
                            coords[dir][targets[dir].back()]);
                }
            }
        }
    }

    PetscFunctionReturn(0);
}

} // end of namespace operators
} // end of namespace petibm
