/***************************************************************************//**
 * \file createR.cpp
 * \author Anush Krishnan (anus@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Definition of function `createR`.
 */


// STL
# include <vector>
# include <functional>

// here goes PETSc headers
# include <petscmat.h>

// here goes headers from our PetIBM
# include "CartesianMesh.h"


/** \brief a short name for the signatures of kernels. 
 * 
 * A kernel is a function that can be used to calculate the diagonal values 
 * according to a given set of i, j, k index
 */
typedef std::function<PetscReal(
        const PetscInt &, const PetscInt &, const PetscInt &)>  KernelType;


/**
 * \brief a function to create diagonal matrix according to input kernel.
 *
 * \param mesh an instance of CartesianMesh.
 * \param kernels a length 3 STL vector holding KernelType.
 * \param D the returned matrix.
 *
 * This is not designed for public use. It's only valid for the functions 
 * defined in this source file.
 *
 * \return PetscErrorCode
 */
PetscErrorCode createDiagMatrix(const CartesianMesh &mesh, 
        const std::vector<KernelType> &kernels, Mat &M)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;
    
    ISLocalToGlobalMapping      *mapping;

    // create matrix
    ierr = DMCreateMatrix(mesh.qPack, &M); CHKERRQ(ierr);
    ierr = MatSeqAIJSetPreallocation(M, 1, nullptr); CHKERRQ(ierr);
    ierr = MatMPIAIJSetPreallocation(M, 1, nullptr, 0, nullptr); CHKERRQ(ierr);
    ierr = MatSetOption(M, MAT_NO_OFF_PROC_ENTRIES, PETSC_TRUE); CHKERRQ(ierr);
    ierr = MatSetOption(M, MAT_KEEP_NONZERO_PATTERN, PETSC_FALSE); CHKERRQ(ierr);
    ierr = MatSetOption(M, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE); CHKERRQ(ierr);

    // get the mapping for each sub-DM
    ierr = DMCompositeGetISLocalToGlobalMappings(mesh.qPack, &mapping); CHKERRQ(ierr);

    // set values to matrix
    for(PetscInt field=0; field<mesh.dim; ++field)
        for(PetscInt k=mesh.bg[field][2]; k<mesh.ed[field][2]; ++k)
            for(PetscInt j=mesh.bg[field][1]; j<mesh.ed[field][1]; ++j)
                for(PetscInt i=mesh.bg[field][0]; i<mesh.ed[field][0]; ++i)
                {
                    PetscInt    idx;
                    MatStencil  self = {k, j, i, 1};

                    // get index local to this process
                    ierr = DMDAConvertToCell(
                            mesh.da[field], self, &idx); CHKERRQ(ierr);

                    // map to global index
                    ierr = ISLocalToGlobalMappingApply(
                            mapping[field], 1, &idx, &idx); CHKERRQ(ierr);

                    // set value
                    ierr = MatSetValue(M, idx, idx, kernels[field](i, j, k), 
                            INSERT_VALUES); CHKERRQ(ierr);
                }

    // assemble matrix
    ierr = MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    // destroy mapping
    for(PetscInt i=0; i<mesh.dim; ++i)
    {
        ierr = ISLocalToGlobalMappingDestroy(&mapping[i]); CHKERRQ(ierr);
    }
    ierr = PetscFree(mapping); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/** \copydoc createR(const CartesianMesh &, Mat &). */
PetscErrorCode createR(const CartesianMesh &mesh, Mat &R)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;
    
    // kernels are kernel functions to calculate entry values
    std::vector<KernelType>     kernel(3);

    // set kernels
    kernel[0] = [&mesh](const PetscInt &i, const PetscInt &j, const PetscInt &k)
        -> PetscReal { return mesh.dL[0][1][j] * mesh.dL[0][2][k]; };
    kernel[1] = [&mesh](const PetscInt &i, const PetscInt &j, const PetscInt &k) 
        -> PetscReal { return mesh.dL[1][0][i] * mesh.dL[1][2][k]; };
    kernel[2] = [&mesh](const PetscInt &i, const PetscInt &j, const PetscInt &k) 
        -> PetscReal { return mesh.dL[2][0][i] * mesh.dL[2][1][j]; };

    ierr = createDiagMatrix(mesh, kernel, R); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/** \copydoc createRInv(const CartesianMesh &, Mat &). */
PetscErrorCode createRInv(const CartesianMesh &mesh, Mat &RInv)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;
    
    // kernels are kernel functions to calculate entry values
    std::vector<KernelType>     kernel(3);

    // set kernels
    kernel[0] = [&mesh](const PetscInt &i, const PetscInt &j, const PetscInt &k)
        -> PetscReal { return 1.0 / (mesh.dL[0][1][j] * mesh.dL[0][2][k]); };
    kernel[1] = [&mesh](const PetscInt &i, const PetscInt &j, const PetscInt &k) 
        -> PetscReal { return 1.0 / (mesh.dL[1][0][i] * mesh.dL[1][2][k]); };
    kernel[2] = [&mesh](const PetscInt &i, const PetscInt &j, const PetscInt &k) 
        -> PetscReal { return 1.0 / (mesh.dL[2][0][i] * mesh.dL[2][1][j]); };

    ierr = createDiagMatrix(mesh, kernel, RInv); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
