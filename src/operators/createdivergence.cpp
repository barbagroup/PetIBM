/**
 * \file createdivergence.cpp
 * \brief Definition of functions creating divergence operator.
 * \author Anush Krishnan (anus@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \copyright MIT.
 */


// STL
# include <functional>

// here goes PETSc headers
# include <petscmat.h>

// here goes headers from our PetIBM
# include <petibm/type.h>
# include <petibm/mesh.h>
# include <petibm/boundary.h>


namespace petibm
{
namespace operators
{

/** \brief a private struct for matrix-free constraint mat. */
struct DivergenceCtx
{
    const type::Boundary    bc;
    type::MatrixModifier    modifier;
    
    DivergenceCtx(const type::Boundary &_bc): 
        bc(_bc), modifier(bc->dim) {};
};


/** \brief a user-defined MatMult for constraint mat. */
PetscErrorCode DCorrectionMult(Mat mat, Vec x, Vec y);


/** \brief a user-defined MatDestroy for constraint mat. */
PetscErrorCode DCorrectionDestroy(Mat mat);


/** \brief a function type for functions returning neighbor's stencil. */
typedef std::function<MatStencil(
        const PetscInt &, const PetscInt &, const PetscInt &)> GetNeighborFunc;


/** \brief a function type for functions returning the value of a matrix entry. */
typedef std::function<PetscReal(
        const PetscInt &, const PetscInt &, const PetscInt &)> Kernel;


// implementation of petibm::operators::createDivergence
PetscErrorCode createDivergence(const type::Mesh &mesh,
                                const type::Boundary &bc, 
                                Mat &D,
                                Mat &DCorrection,
                                const PetscBool &normalize)
{
    PetscFunctionBeginUser;

    PetscErrorCode                  ierr;
    
    std::vector<GetNeighborFunc>    getNeighbor(3);

    std::vector<Kernel>             kernel(3);

    DivergenceCtx                   *ctx;


    // initialize modifier, regardless what's inside it now
    ctx = new DivergenceCtx(bc);


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
            -> PetscReal { return mesh->dL[0][1][j] * mesh->dL[0][2][k]; };
        kernel[1] = [&mesh](const PetscInt &i, const PetscInt &j, const PetscInt &k)
            -> PetscReal { return mesh->dL[1][0][i] * mesh->dL[1][2][k]; };
        kernel[2] = [&mesh](const PetscInt &i, const PetscInt &j, const PetscInt &k)
            -> PetscReal { return mesh->dL[2][0][i] * mesh->dL[2][1][j]; };
    }


    // create matrix
    ierr = MatCreate(mesh->comm, &D); CHKERRQ(ierr);
    ierr = MatSetSizes(D, mesh->pNLocal, mesh->UNLocal, 
            PETSC_DETERMINE, PETSC_DETERMINE); CHKERRQ(ierr);
    ierr = MatSetFromOptions(D); CHKERRQ(ierr);
    ierr = MatSeqAIJSetPreallocation(D, mesh->dim*2, nullptr); CHKERRQ(ierr);

    // TODO: a better guess for allocation in MPIAIJ
    ierr = MatMPIAIJSetPreallocation(D, mesh->dim*2, nullptr, 3, nullptr); CHKERRQ(ierr);

    ierr = MatSetUp(D); CHKERRQ(ierr);
    ierr = MatSetOption(D, MAT_KEEP_NONZERO_PATTERN, PETSC_FALSE); CHKERRQ(ierr);
    ierr = MatSetOption(D, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE); CHKERRQ(ierr);


    // set values to matrix
    for(PetscInt k=mesh->bg[3][2]; k<mesh->ed[3][2]; ++k)
        for(PetscInt j=mesh->bg[3][1]; j<mesh->ed[3][1]; ++j)
            for(PetscInt i=mesh->bg[3][0]; i<mesh->ed[3][0]; ++i)
            {
                PetscInt    self;
                
                // row index is based on pressure grid
                ierr = mesh->getGlobalIndex(3, {k, j, i, 0}, self); CHKERRQ(ierr);

                // calculate and set values of the two surfaces in each direction
                for(PetscInt field=0; field<mesh->dim; ++field)
                {
                    // the stencil of target surface
                    MatStencil      target;

                    // index of target surface in packed Vec
                    PetscInt        targetIdx;

                    // the absolute value of coefficient
                    PetscReal       value = kernel[field](i, j, k);


                    // the surface toward positive direction
                    target = {k, j, i, 0};

                    ierr = mesh->getPackedGlobalIndex(field, target, targetIdx); CHKERRQ(ierr);

                    ierr = MatSetValue(D, self, targetIdx, 
                            value, INSERT_VALUES); CHKERRQ(ierr);

                    // store the coefficient of ghost points for future use
                    if (targetIdx == -1) 
                        ctx->modifier[field][target] = {self, value};


                    // the surface toward negative direction
                    target = getNeighbor[field](i, j, k);

                    ierr = mesh->getPackedGlobalIndex(field, target, targetIdx); CHKERRQ(ierr);

                    // note the value for this face is just a negative version
                    ierr = MatSetValue(D, self, targetIdx, 
                            - value, INSERT_VALUES); CHKERRQ(ierr);

                    // store the coefficient of ghost points for future use
                    if (targetIdx == -1) 
                        ctx->modifier[field][target] = {self, - value};
                }
            }


    // temporarily assemble matrix
    ierr = MatAssemblyBegin(D, MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(D, MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);


    // modify the column of the neighbors of ghost points
    // Only ghost points on Neumann BC will modify Divergence operator
    for(PetscInt f=0; f<mesh->dim; ++f)
        for(auto &bd: bc->bds[f])
            if (bd->onThisProc)
                for(auto &pt: bd->points)
                {
                    PetscInt    col = pt.second.targetPackedId;
                    PetscReal   value =
                        ctx->modifier[f][pt.first].coeff * pt.second.a0;

                    ierr = MatSetValue(D, ctx->modifier[f][pt.first].row,
                            col, value, ADD_VALUES); CHKERRQ(ierr);
                }


    // assemble matrix, an implicit mpi barrier is applied
    ierr = MatAssemblyBegin(D, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(D, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);


    // create a matrix-free constraint matrix for boundary correction
    ierr = MatCreateShell(mesh->comm, mesh->pNLocal, mesh->UNLocal,
            PETSC_DETERMINE, PETSC_DETERMINE, (void *) ctx, &DCorrection); 
    CHKERRQ(ierr);

    ierr = MatShellSetOperation(DCorrection, MATOP_MULT,
            (void(*)(void)) DCorrectionMult); CHKERRQ(ierr);

    ierr = MatShellSetOperation(DCorrection, MATOP_DESTROY,
            (void(*)(void)) DCorrectionDestroy); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


// implementation of DCorrectionMult
PetscErrorCode DCorrectionMult(Mat mat, Vec x, Vec y)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    DivergenceCtx       *ctx;

    // get the context
    ierr = MatShellGetContext(mat, (void *) &ctx); CHKERRQ(ierr);

    // zero the output vector
    ierr = VecSet(y, 0.0); CHKERRQ(ierr);

    // set the correction values to corresponding rows
    for(PetscInt f=0; f<ctx->bc->dim; ++f)
        for(auto &bd: ctx->bc->bds[f])
            if (bd->onThisProc)
                for(auto &pt: bd->points)
                {
                    ierr = VecSetValue(y, ctx->modifier[f][pt.first].row,
                            ctx->modifier[f][pt.first].coeff * pt.second.a1, 
                            ADD_VALUES); CHKERRQ(ierr);
                }

    // assembly (not sure if this is necessary and if this causes overhead)
    // but there is an implicit MPI barrier in assembly, which we defenitely need
    ierr = VecAssemblyBegin(y); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(y); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


// implementation of DCorrectionDestroy
PetscErrorCode DCorrectionDestroy(Mat mat)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;
    DivergenceCtx       *ctx;

    // get the context
    ierr = MatShellGetContext(mat, (void *) &ctx); CHKERRQ(ierr);

    // deallocate the context
    delete ctx;

    PetscFunctionReturn(0);
}

} // end of namespace operators
} // end of namespace petibm
