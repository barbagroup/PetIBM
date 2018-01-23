/**
 * \file creatediagmatrix.cpp
 * \brief Definitions of functions creating different kinds of diagonal matrices.
 * \author Anush Krishnan (anush@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \copyright MIT.
 */


// STL
# include <vector>
# include <functional>

// here goes PETSc headers
# include <petscmat.h>

// here goes headers from our PetIBM
# include <petibm/mesh.h>


namespace petibm
{
namespace operators
{

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
 * \param M the returned matrix.
 *
 * This is not designed for public use. It's only valid for the functions 
 * defined in this source file.
 *
 * \return PetscErrorCode
 */
PetscErrorCode createDiagMatrix(const type::Mesh &mesh, 
                                const std::vector<KernelType> &kernels,
                                Mat &M)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;
    
    // create matrix
    ierr = MatCreate(mesh->comm, &M); CHKERRQ(ierr);
    ierr = MatSetSizes(M, mesh->UNLocal, mesh->UNLocal, 
            PETSC_DETERMINE, PETSC_DETERMINE); CHKERRQ(ierr);
    ierr = MatSetFromOptions(M); CHKERRQ(ierr);
    ierr = MatSeqAIJSetPreallocation(M, 1, nullptr); CHKERRQ(ierr);
    ierr = MatMPIAIJSetPreallocation(M, 1, nullptr, 0, nullptr); CHKERRQ(ierr);
    ierr = MatSetUp(M); CHKERRQ(ierr);
    ierr = MatSetOption(M, MAT_NO_OFF_PROC_ENTRIES, PETSC_TRUE); CHKERRQ(ierr);
    ierr = MatSetOption(M, MAT_KEEP_NONZERO_PATTERN, PETSC_FALSE); CHKERRQ(ierr);
    ierr = MatSetOption(M, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE); CHKERRQ(ierr);


    // set values to matrix
    for(PetscInt field=0; field<mesh->dim; ++field)
        for(PetscInt k=mesh->bg[field][2]; k<mesh->ed[field][2]; ++k)
            for(PetscInt j=mesh->bg[field][1]; j<mesh->ed[field][1]; ++j)
                for(PetscInt i=mesh->bg[field][0]; i<mesh->ed[field][0]; ++i)
                {
                    PetscInt    idx;
                    MatStencil  self = {k, j, i, 1};

                    // get packed index of this velocity point
                    ierr = mesh->getPackedGlobalIndex(field, self, idx); CHKERRQ(ierr);

                    // set value
                    ierr = MatSetValue(M, idx, idx, kernels[field](i, j, k), 
                            INSERT_VALUES); CHKERRQ(ierr);
                }

    // assemble matrix
    ierr = MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // createDiagMatrix


// implementation of petibm::operators::createR
PetscErrorCode createR(const type::Mesh &mesh, Mat &R)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;
    
    // kernels are kernel functions to calculate entry values
    std::vector<KernelType>     kernel(3);

    // set kernels
    kernel[0] = [&mesh](const PetscInt &i, const PetscInt &j, const PetscInt &k)
        -> PetscReal { return mesh->dL[0][1][j] * mesh->dL[0][2][k]; };
    kernel[1] = [&mesh](const PetscInt &i, const PetscInt &j, const PetscInt &k) 
        -> PetscReal { return mesh->dL[1][0][i] * mesh->dL[1][2][k]; };
    kernel[2] = [&mesh](const PetscInt &i, const PetscInt &j, const PetscInt &k) 
        -> PetscReal { return mesh->dL[2][0][i] * mesh->dL[2][1][j]; };

    // call the function to create matrix
    ierr = createDiagMatrix(mesh, kernel, R); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // createR


// implementation of petibm::operators::createRInv
PetscErrorCode createRInv(const type::Mesh &mesh, Mat &RInv)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;
    
    // kernels are kernel functions to calculate entry values
    std::vector<KernelType>     kernel(3);

    // set kernels
    kernel[0] = [&mesh](const PetscInt &i, const PetscInt &j, const PetscInt &k)
        -> PetscReal { return 1.0 / (mesh->dL[0][1][j] * mesh->dL[0][2][k]); };
    kernel[1] = [&mesh](const PetscInt &i, const PetscInt &j, const PetscInt &k) 
        -> PetscReal { return 1.0 / (mesh->dL[1][0][i] * mesh->dL[1][2][k]); };
    kernel[2] = [&mesh](const PetscInt &i, const PetscInt &j, const PetscInt &k) 
        -> PetscReal { return 1.0 / (mesh->dL[2][0][i] * mesh->dL[2][1][j]); };

    // call the function to create matrix
    ierr = createDiagMatrix(mesh, kernel, RInv); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // createRInv


// implementation of petibm::operators::createMHead
PetscErrorCode createMHead(const type::Mesh &mesh, Mat &MHead)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;
    
    // kernels are kernel functions to calculate entry values
    std::vector<KernelType>     kernel(3);

    // set kernels
    kernel[0] = [&mesh](const PetscInt &i, const PetscInt &j, const PetscInt &k)
        -> PetscReal { return mesh->dL[0][0][i]; };
    kernel[1] = [&mesh](const PetscInt &i, const PetscInt &j, const PetscInt &k) 
        -> PetscReal { return mesh->dL[1][1][j]; };
    kernel[2] = [&mesh](const PetscInt &i, const PetscInt &j, const PetscInt &k) 
        -> PetscReal { return mesh->dL[2][2][k]; };

    // call the function to create matrix
    ierr = createDiagMatrix(mesh, kernel, MHead); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // createMHead


// implementation of petibm::operators::createM
PetscErrorCode createM(const type::Mesh &mesh, Mat &M)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;
    
    // kernels are kernel functions to calculate entry values
    std::vector<KernelType>     kernel(3);

    // set kernels
    kernel[0] = [&mesh](const PetscInt &i, const PetscInt &j, const PetscInt &k)
        -> PetscReal { return mesh->dL[0][0][i] / (mesh->dL[0][1][j] * mesh->dL[0][2][k]); };
    kernel[1] = [&mesh](const PetscInt &i, const PetscInt &j, const PetscInt &k) 
        -> PetscReal { return mesh->dL[1][1][j] / (mesh->dL[1][0][i] * mesh->dL[1][2][k]); };
    kernel[2] = [&mesh](const PetscInt &i, const PetscInt &j, const PetscInt &k) 
        -> PetscReal { return mesh->dL[2][2][k] / (mesh->dL[2][0][i] * mesh->dL[2][1][j]); };

    // call the function to create matrix
    ierr = createDiagMatrix(mesh, kernel, M); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // createM


// implementation of petibm::operators::createIdentity
PetscErrorCode createIdentity(const type::Mesh &mesh, Mat &I)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;
    
    // kernels are kernel functions to calculate entry values
    std::vector<KernelType>     kernel(3);

    // set kernels
    kernel[0] = [&mesh](const PetscInt &i, const PetscInt &j, const PetscInt &k)
        -> PetscReal { return 1.0; };
    kernel[1] = [&mesh](const PetscInt &i, const PetscInt &j, const PetscInt &k) 
        -> PetscReal { return 1.0; };
    kernel[2] = [&mesh](const PetscInt &i, const PetscInt &j, const PetscInt &k) 
        -> PetscReal { return 1.0; };

    // call the function to create matrix
    ierr = createDiagMatrix(mesh, kernel, I); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // createIdentity

} // end of namespace operators
} // end of namespace petibm
