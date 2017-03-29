/***************************************************************************//**
 * \file createDivergence.cpp
 * \author Anush Krishnan (anus@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Definition of functions regarding to creating divergence operators.
 */


// here goes PETSc headers
# include <petscmat.h>

// here goes headers from our PetIBM
# include "CartesianMesh.h"
# include "types.h"


/** \brief a function type for functions returning neighbor's stencil. */
typedef std::function<MatStencil(
        const PetscInt &, const PetscInt &, const PetscInt &)> GetNeighborFunc;


/** \brief a function type for functions returning the value of a matrix entry. */
typedef std::function<PetscReal(
        const PetscInt &, const PetscInt &, const PetscInt &)> Kernel;


/** \copydoc createDivergence. */
PetscErrorCode createDivergence(
        const CartesianMesh &mesh, Mat &D, const PetscBool &normalize)
{
    PetscFunctionBeginUser;

    PetscErrorCode                  ierr;
    
    std::vector<GetNeighborFunc>    getNeighbor(3);

    std::vector<Kernel>             kernel(3);


    // set up getNeighbor
    getNeighbor[0] = [](const PetscInt &i, const PetscInt &j, const PetscInt &k)
        -> MatStencil { return {k, j, i-1}; };
    getNeighbor[1] = [](const PetscInt &i, const PetscInt &j, const PetscInt &k)
        -> MatStencil { return {k, j-1, i}; };
    getNeighbor[2] = [](const PetscInt &i, const PetscInt &j, const PetscInt &k)
        -> MatStencil { return {k-1, j, i}; };


    // set up kernel
    if (normalize)
        kernel[0] = kernel[1] = kernel[2] = 
            std::bind([]() -> PetscReal { return 1.0; });
    else
    {
        kernel[0] = [&mesh](const PetscInt &i, const PetscInt &j, const PetscInt &k)
            -> PetscReal { return mesh.dL[0][1][j] * mesh.dL[0][2][k]; };
        kernel[1] = [&mesh](const PetscInt &i, const PetscInt &j, const PetscInt &k)
            -> PetscReal { return mesh.dL[1][0][i] * mesh.dL[1][2][k]; };
        kernel[2] = [&mesh](const PetscInt &i, const PetscInt &j, const PetscInt &k)
            -> PetscReal { return mesh.dL[2][0][i] * mesh.dL[2][1][j]; };
    }


    // create matrix
    ierr = MatCreate(*mesh.comm, &D); CHKERRQ(ierr);
    ierr = MatSetSizes(D, mesh.lambdaNLocal, mesh.qNLocal, 
            PETSC_DETERMINE, PETSC_DETERMINE); CHKERRQ(ierr);
    ierr = MatSetFromOptions(D); CHKERRQ(ierr);
    ierr = MatSeqAIJSetPreallocation(D, mesh.dim*2, nullptr); CHKERRQ(ierr);

    // TODO: a better guess for allocation in MPIAIJ
    ierr = MatMPIAIJSetPreallocation(D, mesh.dim*2, nullptr, 3, nullptr); CHKERRQ(ierr);

    ierr = MatSetUp(D); CHKERRQ(ierr);
    ierr = MatSetOption(D, MAT_KEEP_NONZERO_PATTERN, PETSC_FALSE); CHKERRQ(ierr);
    ierr = MatSetOption(D, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE); CHKERRQ(ierr);


    // set values to matrix
    for(PetscInt k=mesh.bg[3][2]; k<mesh.ed[3][2]; ++k)
        for(PetscInt j=mesh.bg[3][1]; j<mesh.ed[3][1]; ++j)
            for(PetscInt i=mesh.bg[3][0]; i<mesh.ed[3][0]; ++i)
            {
                PetscInt    idx, self;

                ierr = DMDAConvertToCell(mesh.da[3], {k, j, i, 0}, &self);
                CHKERRQ(ierr);

                ierr = ISLocalToGlobalMappingApply(
                        mesh.lambdaMapping[0], 1, &self, &self);
                CHKERRQ(ierr);

                for(PetscInt field=0; field<mesh.dim; ++field)
                {
                    ierr = DMDAConvertToCell(
                            mesh.da[field], {k, j, i, 0}, &idx); CHKERRQ(ierr);

                    ierr = ISLocalToGlobalMappingApply(
                            mesh.qMapping[field], 1, &idx, &idx); CHKERRQ(ierr);

                    ierr = MatSetValue(D, self, idx, 
                            kernel[field](i, j, k), INSERT_VALUES); CHKERRQ(ierr);

                    ierr = DMDAConvertToCell(mesh.da[field], 
                            getNeighbor[field](i, j, k), &idx); CHKERRQ(ierr);

                    ierr = ISLocalToGlobalMappingApply(
                            mesh.qMapping[field], 1, &idx, &idx); CHKERRQ(ierr);

                    ierr = MatSetValue(D, self, idx, 
                            - kernel[field](i, j, k), INSERT_VALUES); CHKERRQ(ierr);
                }
            }


    // TODO: modify diag values for Neumann BCs.


    // assemble matrix
    ierr = MatAssemblyBegin(D, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(D, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
