/***************************************************************************//**
 * \file createGradient.cpp
 * \author Anush Krishnan (anus@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Definition of functions regarding to creating gradient operators.
 */


// here goes PETSc headers
# include <petscmat.h>

// here goes headers from our PetIBM
# include "CartesianMesh.h"
# include "types.h"


/** \brief STL vector holding MatStencil. */
typedef std::vector<MatStencil> StencilVec;


/** \brief a function type for functions returning neighbors' stencils. */
typedef std::function<StencilVec(
        const PetscInt &, const PetscInt &, const PetscInt &)> GetNeighborFunc;


/** \brief a function type for functions computing matrix entries' values. */
typedef std::function<types::RealVec1D(
        const PetscInt &, const PetscInt &, const PetscInt &)> Kernel;


/** \copydoc createGradient. */
PetscErrorCode createGradient(
        const CartesianMesh &mesh, Mat &G, const PetscBool &normalize)
{
    PetscFunctionBeginUser;

    PetscErrorCode                  ierr;
    
    ISLocalToGlobalMapping          *qMapping, *pMapping;

    std::vector<GetNeighborFunc>    getNeighbor(3);

    std::vector<Kernel>             kernel(3);


    // total number of local velocity points
    PetscInt                nQLcl =
        mesh.m[0][0] * mesh.m[0][1] * mesh.m[0][2] + 
        mesh.m[1][0] * mesh.m[1][1] * mesh.m[1][2] + 
        mesh.m[2][0] * mesh.m[2][1] * mesh.m[2][2];

    // total number of local pressure points
    PetscInt                nPLcl = mesh.m[3][0] * mesh.m[3][1] * mesh.m[3][2];


    // set up the function used to get the IDs of neighbors
    getNeighbor[0] = [](const PetscInt &i, const PetscInt &j, const PetscInt &k)
        -> std::vector<MatStencil> { return {{k, j, i, 0}, {k, j, i+1, 0}}; };

    getNeighbor[1] = [](const PetscInt &i, const PetscInt &j, const PetscInt &k)
        -> std::vector<MatStencil> { return {{k, j, i, 0}, {k, j+1, i, 0}}; };

    getNeighbor[2] = [](const PetscInt &i, const PetscInt &j, const PetscInt &k)
        -> std::vector<MatStencil> { return {{k, j, i, 0}, {k+1, j, i, 0}}; };


    // set up the kernel that calculates the entry values
    if (normalize)
        kernel[0] = kernel[1] = kernel[2] = std::bind(
                []() -> types::RealVec1D { return {-1.0, 1.0}; });
    else
    {
        kernel[0] = [&mesh](const PetscInt &i, const PetscInt &j, const PetscInt &k)
            -> types::RealVec1D {PetscReal v=1.0/mesh.dL[0][0][i]; return {-v, v};};

        kernel[1] = [&mesh](const PetscInt &i, const PetscInt &j, const PetscInt &k)
            -> types::RealVec1D {PetscReal v=1.0/mesh.dL[1][1][j]; return {-v, v};};

        kernel[2] = [&mesh](const PetscInt &i, const PetscInt &j, const PetscInt &k)
            -> types::RealVec1D {PetscReal v=1.0/mesh.dL[2][2][k]; return {-v, v};};
    }


    // create matrix
    ierr = MatCreate(*mesh.comm, &G); CHKERRQ(ierr);
    ierr = MatSetSizes(
            G, nQLcl, nPLcl, PETSC_DETERMINE, PETSC_DETERMINE); CHKERRQ(ierr);
    ierr = MatSetFromOptions(G); CHKERRQ(ierr);
    ierr = MatSeqAIJSetPreallocation(G, 2, nullptr); CHKERRQ(ierr);
    ierr = MatMPIAIJSetPreallocation(G, 2, nullptr, 1, nullptr); CHKERRQ(ierr);
    ierr = MatSetUp(G); CHKERRQ(ierr);
    ierr = MatSetOption(G, MAT_KEEP_NONZERO_PATTERN, PETSC_FALSE); CHKERRQ(ierr);
    ierr = MatSetOption(G, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE); CHKERRQ(ierr);


    // get the mapping for each sub-DM in qPack
    ierr = DMCompositeGetISLocalToGlobalMappings(
            mesh.qPack, &qMapping); CHKERRQ(ierr);

    // get the mapping for each sub-DM in lambdaPack
    ierr = DMCompositeGetISLocalToGlobalMappings(
            mesh.lambdaPack, &pMapping); CHKERRQ(ierr);


    // set values to matrix
    for(PetscInt field=0; field<mesh.dim; ++field)
        for(PetscInt k=mesh.bg[field][2]; k<mesh.ed[field][2]; ++k)
            for(PetscInt j=mesh.bg[field][1]; j<mesh.ed[field][1]; ++j)
                for(PetscInt i=mesh.bg[field][0]; i<mesh.ed[field][0]; ++i)
                {
                    PetscInt            rId, cId;

                    StencilVec          loc = getNeighbor[field](i, j, k);

                    types::RealVec1D    values = kernel[field](i, j, k);

                    // get velocity index local to this process
                    ierr = DMDAConvertToCell(mesh.da[field], loc[0], &rId); 
                    CHKERRQ(ierr);

                    // map local velocity index to global index
                    ierr = ISLocalToGlobalMappingApply(
                            qMapping[field], 1, &rId, &rId); CHKERRQ(ierr);

                    // loop through columns
                    for(PetscInt n=0; n<2; ++n)
                    {
                        ierr = DMDAConvertToCell(
                                mesh.da[3], loc[n], &cId); CHKERRQ(ierr);

                        ierr = ISLocalToGlobalMappingApply(
                                pMapping[0], 1, &cId, &cId); CHKERRQ(ierr);

                        ierr = MatSetValue(G, rId, cId, 
                                values[n], INSERT_VALUES); CHKERRQ(ierr);
                    }
                }


    // assemble matrix
    ierr = MatAssemblyBegin(G, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(G, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);


    // destroy qMapping
    for(PetscInt i=0; i<mesh.dim; ++i)
    {
        ierr = ISLocalToGlobalMappingDestroy(&qMapping[i]); CHKERRQ(ierr);
    }
    ierr = PetscFree(qMapping); CHKERRQ(ierr);

    // get the number of sub-DMs in lambdaPack for the following loop
    ierr = DMCompositeGetNumberDM(mesh.lambdaPack, &nPLcl); CHKERRQ(ierr);

    // destroy pMapping
    for(PetscInt i=0; i<nPLcl; ++i)
    {
        ierr = ISLocalToGlobalMappingDestroy(&pMapping[i]); CHKERRQ(ierr);
    }
    ierr = PetscFree(pMapping); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
