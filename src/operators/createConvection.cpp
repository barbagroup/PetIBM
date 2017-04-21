/***************************************************************************//**
 * \file createConvection.cpp
 * \author Anush Krishnan (anus@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Definition of functions regarding to creating Convection operators.
 */


// TODO: investigate if exact interpolation is necessary


// STL
# include <memory>

// PETSc
# include <petscmat.h>

// PetIBM
# include "CartesianMesh.h"
# include "Boundary.h"


/** \brief a private struct used in MatShell. */
struct NonLinearCtx
{
    std::shared_ptr<const CartesianMesh>    mesh;
    std::shared_ptr<const Boundary>         bc;
    std::vector<Vec>                        qLocal;
};

/**
 * \brief a private function for convection operator's MatMult in 2D.
 *
 * For user-defined PETSc Mat operations, please refer to PETSc manual.
 */
PetscErrorCode ConvectionMult2D(Mat mat, Vec x, Vec y);

/**
 * \brief a private function for convection operator's MatMult in 3D.
 *
 * For user-defined PETSc Mat operations, please refer to PETSc manual.
 */
PetscErrorCode ConvectionMult3D(Mat mat, Vec x, Vec y);

/** \brief user-defined destroyer for convection operator. */
PetscErrorCode ConvectionDestroy(Mat mat);


/** \brief a private kernel for the convection at a u-velocity point in 2D. */
inline PetscReal kernelU(
        NonLinearCtx const * const &ctx, const std::vector<PetscReal**> &flux, 
        const PetscInt &i, const PetscInt &j);


/** \brief a private kernel for the convection at a u-velocity point in 3D. */
inline PetscReal kernelU(
        NonLinearCtx const * const &ctx, const std::vector<PetscReal***> &flux, 
        const PetscInt &i, const PetscInt &j, const PetscInt &k);


/** \brief a private kernel for the convection at a v-velocity point in 2D. */
inline PetscReal kernelV(
        NonLinearCtx const * const &ctx, const std::vector<PetscReal**> &flux, 
        const PetscInt &i, const PetscInt &j);


/** \brief a private kernel for the convection at a v-velocity point in 3D. */
inline PetscReal kernelV(
        NonLinearCtx const * const &ctx, const std::vector<PetscReal***> &flux, 
        const PetscInt &i, const PetscInt &j, const PetscInt &k);


/** \brief a private kernel for the convection at a w-velocity point in 3D. */
inline PetscReal kernelW(
        NonLinearCtx const * const &ctx, const std::vector<PetscReal***> &flux, 
        const PetscInt &i, const PetscInt &j, const PetscInt &k);


/** \copydoc createConvection. */
PetscErrorCode createConvection(
        const CartesianMesh &mesh, const Boundary &bd, Mat &H)
{
    PetscFunctionBeginUser;

    PetscErrorCode  ierr;

    NonLinearCtx    *ctx;

    // global dimension
    PetscInt    N = 
        mesh.n[0][0] * mesh.n[0][1] * mesh.n[0][2] + 
        mesh.n[1][0] * mesh.n[1][1] * mesh.n[1][2] + 
        ((mesh.dim==3)? mesh.n[2][0] * mesh.n[2][1] * mesh.n[2][2] : 0);

    // allocate space for ctx
    ctx = new NonLinearCtx;
    ctx->mesh = std::shared_ptr<const CartesianMesh>(&mesh, [](const CartesianMesh*){});
    ctx->bc = std::shared_ptr<const Boundary>(&bd, [](const Boundary*){});
    ctx->qLocal.resize(ctx->mesh->dim);
    
    // create necessary local vectors
    for(PetscInt f=0; f<ctx->mesh->dim; ++f)
    {
        ierr = DMCreateLocalVector(ctx->mesh->da[f], &ctx->qLocal[f]);
        CHKERRQ(ierr);
    }

    // create a matrix-free operator
    ierr = MatCreateShell(*mesh.comm, mesh.UNLocal, mesh.UNLocal, 
            N, N, (void *) ctx, &H); CHKERRQ(ierr);

    // bind MatMult
    if (ctx->mesh->dim == 2)
    {
        ierr = MatShellSetOperation(H, MATOP_MULT, 
                (void(*)(void)) ConvectionMult2D); CHKERRQ(ierr);
    }
    else // assuem the dim is either 2 or 3
    {
        ierr = MatShellSetOperation(H, MATOP_MULT, 
                (void(*)(void)) ConvectionMult3D); CHKERRQ(ierr);
    }

    // bind MatDestroy
    ierr = MatShellSetOperation(H, MATOP_DESTROY, 
            (void(*)(void)) ConvectionDestroy); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/** \copydoc UserMult2D. */
PetscErrorCode ConvectionMult2D(Mat mat, Vec x, Vec y)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    NonLinearCtx        *ctx;

    std::vector<Vec>    unPacked(2);

    std::vector<PetscReal**>    xArry(2);

    std::vector<PetscReal**>    yArry(2);


    // get the context
    ierr = MatShellGetContext(mat, (void *) &ctx); CHKERRQ(ierr);

    // get local (including overlapped points) values of x
    ierr = DMCompositeScatterArray(ctx->mesh->UPack, 
            x, ctx->qLocal.data()); CHKERRQ(ierr);

    // set the values of ghost points in local vectors
    ierr = ctx->bc->copyValues2LocalVecs(ctx->qLocal); CHKERRQ(ierr);

    // get unPacked vectors of y
    ierr = DMCompositeGetAccessArray(ctx->mesh->UPack,
            y, ctx->mesh->dim, nullptr, unPacked.data()); CHKERRQ(ierr);


    // get underlying data of Vecs
    for(PetscInt f=0; f<ctx->mesh->dim; ++f)
    {
        ierr = DMDAVecGetArrayRead(
                ctx->mesh->da[f], ctx->qLocal[f], &xArry[f]); CHKERRQ(ierr);

        ierr = DMDAVecGetArray(
                ctx->mesh->da[f], unPacked[f], &yArry[f]); CHKERRQ(ierr);
    }


    // u-velocity field
    for(PetscInt j=ctx->mesh->bg[0][1]; j<ctx->mesh->ed[0][1]; ++j)
        for(PetscInt i=ctx->mesh->bg[0][0]; i<ctx->mesh->ed[0][0]; ++i)
            yArry[0][j][i] = kernelU(ctx, xArry, i, j);


    // v-velocity field
    for(PetscInt j=ctx->mesh->bg[1][1]; j<ctx->mesh->ed[1][1]; ++j)
        for(PetscInt i=ctx->mesh->bg[1][0]; i<ctx->mesh->ed[1][0]; ++i)
            yArry[1][j][i] = kernelV(ctx, xArry, i, j);


    // return underlying arrays of Vecs
    for(PetscInt f=0; f<ctx->mesh->dim; ++f)
    {
        ierr = DMDAVecRestoreArrayRead(
                ctx->mesh->da[f], ctx->qLocal[f], &xArry[f]); CHKERRQ(ierr);
    
        ierr = DMDAVecRestoreArray(
                ctx->mesh->da[f], unPacked[f], &yArry[f]); CHKERRQ(ierr);
    }


    // return the ownership of unpacked vectors back to packed vectors
    ierr = DMCompositeRestoreAccessArray(ctx->mesh->UPack,
            y, ctx->mesh->dim, nullptr, unPacked.data()); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/** \copydoc UserMult3D. */
PetscErrorCode ConvectionMult3D(Mat mat, Vec x, Vec y)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    NonLinearCtx        *ctx;

    std::vector<Vec>    unPacked(3);

    std::vector<PetscReal***>    xArry(3);

    std::vector<PetscReal***>    yArry(3);


    // get the context
    ierr = MatShellGetContext(mat, (void *) &ctx); CHKERRQ(ierr);

    // get local (including overlapped points) vectors of x
    ierr = DMCompositeScatterArray(ctx->mesh->UPack, 
            x, ctx->qLocal.data()); CHKERRQ(ierr);

    // set the values of ghost points in local vectors
    ierr = ctx->bc->copyValues2LocalVecs(ctx->qLocal); CHKERRQ(ierr);

    // get unPacked vectors of y
    ierr = DMCompositeGetAccessArray(ctx->mesh->UPack,
            y, ctx->mesh->dim, nullptr, unPacked.data()); CHKERRQ(ierr);


    // get underlying data of Vecs
    for(PetscInt f=0; f<ctx->mesh->dim; ++f)
    {
        ierr = DMDAVecGetArrayRead(
                ctx->mesh->da[f], ctx->qLocal[f], &xArry[f]); CHKERRQ(ierr);

        ierr = DMDAVecGetArray(
                ctx->mesh->da[f], unPacked[f], &yArry[f]); CHKERRQ(ierr);
    }


    // u-velocity field
    for(PetscInt k=ctx->mesh->bg[0][2]; k<ctx->mesh->ed[0][2]; ++k)
        for(PetscInt j=ctx->mesh->bg[0][1]; j<ctx->mesh->ed[0][1]; ++j)
            for(PetscInt i=ctx->mesh->bg[0][0]; i<ctx->mesh->ed[0][0]; ++i)
                yArry[0][k][j][i] = kernelU(ctx, xArry, i, j, k);


    // v-velocity field
    for(PetscInt k=ctx->mesh->bg[1][2]; k<ctx->mesh->ed[1][2]; ++k)
        for(PetscInt j=ctx->mesh->bg[1][1]; j<ctx->mesh->ed[1][1]; ++j)
            for(PetscInt i=ctx->mesh->bg[1][0]; i<ctx->mesh->ed[1][0]; ++i)
                yArry[1][k][j][i] = kernelV(ctx, xArry, i, j, k);


    // w-velocity field
    for(PetscInt k=ctx->mesh->bg[2][2]; k<ctx->mesh->ed[2][2]; ++k)
        for(PetscInt j=ctx->mesh->bg[2][1]; j<ctx->mesh->ed[2][1]; ++j)
            for(PetscInt i=ctx->mesh->bg[2][0]; i<ctx->mesh->ed[2][0]; ++i)
                yArry[2][k][j][i] = kernelW(ctx, xArry, i, j, k);


    // return underlying arrays of Vecs
    for(PetscInt f=0; f<ctx->mesh->dim; ++f)
    {
        ierr = DMDAVecRestoreArrayRead(
                ctx->mesh->da[f], ctx->qLocal[f], &xArry[f]); CHKERRQ(ierr);
    
        ierr = DMDAVecRestoreArray(
                ctx->mesh->da[f], unPacked[f], &yArry[f]); CHKERRQ(ierr);
    }


    // return the ownership of unpacked vectors back to packed vectors
    ierr = DMCompositeRestoreAccessArray(ctx->mesh->UPack,
            y, ctx->mesh->dim, nullptr, unPacked.data()); CHKERRQ(ierr);


    PetscFunctionReturn(0);
}


/** \copydoc UserDestroy. */
PetscErrorCode ConvectionDestroy(Mat mat)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    NonLinearCtx        *ctx;

    // get the context
    ierr = MatShellGetContext(mat, (void *) &ctx); CHKERRQ(ierr);

    // destroy qLocal
    for(PetscInt f=0; f<ctx->mesh->dim; f++)
    {
        ierr = VecDestroy(&ctx->qLocal[f]); CHKERRQ(ierr);
    }

    // deallocate the memory space pointed by ctx
    delete ctx;

    PetscFunctionReturn(0);
}


/** \brief a private kernel for the convection at a u-velocity point in 2D. */
inline PetscReal kernelU(
        NonLinearCtx const * const &ctx, const std::vector<PetscReal**> &flux, 
        const PetscInt &i, const PetscInt &j)
{
    PetscReal   uSelf;
    PetscReal   uS, uN, uW, uE;
    PetscReal   vS, vN;

    // prepare self
    uSelf = flux[0][j][i];

    // prepare u
    uW = (uSelf + flux[0][j][i-1]) / 2.0;
    uE = (uSelf + flux[0][j][i+1]) / 2.0;
    uS = (uSelf + flux[0][j-1][i]) / 2.0;
    uN = (uSelf + flux[0][j+1][i]) / 2.0;

    // prepare v
    vS = (flux[1][j-1][i] + flux[1][j-1][i+1]) / 2.0;
    vN = (flux[1][j][i] + flux[1][j][i+1]) / 2.0;

    return 
        (uE * uE - uW * uW) / ctx->mesh->dL[0][0][i] + 
        (vN * uN - vS * uS) / ctx->mesh->dL[0][1][j];
}


/** \brief a private kernel for the convection at a v-velocity point in 2D. */
inline PetscReal kernelV(
        NonLinearCtx const * const &ctx, const std::vector<PetscReal**> &flux, 
        const PetscInt &i, const PetscInt &j)
{
    PetscReal   vSelf;
    PetscReal   uW, uE;
    PetscReal   vS, vN, vW, vE;

    // prepare self
    vSelf = flux[1][j][i];

    // prepare u
    uW = (flux[0][j][i-1] + flux[0][j+1][i-1]) / 2.0;
    uE = (flux[0][j][i] + flux[0][j+1][i]) / 2.0;

    // prepare v
    vW = (vSelf + flux[1][j][i-1]) / 2.0;
    vE = (vSelf + flux[1][j][i+1]) / 2.0;
    vS = (vSelf + flux[1][j-1][i]) / 2.0;
    vN = (vSelf + flux[1][j+1][i]) / 2.0;

    return 
        (uE * vE - uW * vW) / ctx->mesh->dL[1][0][i] +
        (vN * vN - vS * vS) / ctx->mesh->dL[1][1][j];
}


/** \brief a private kernel for the convection at a u-velocity point in 3D. */
inline PetscReal kernelU(
        NonLinearCtx const * const &ctx, const std::vector<PetscReal***> &flux, 
        const PetscInt &i, const PetscInt &j, const PetscInt &k)
{
    PetscReal   uSelf;
    PetscReal   uS, uN, uW, uE, uB, uF;
    PetscReal   vS, vN;
    PetscReal   wB, wF;

    // prepare self
    uSelf = flux[0][k][j][i];

    // prepare u
    uW = (uSelf + flux[0][k][j][i-1]) / 2.0;
    uE = (uSelf + flux[0][k][j][i+1]) / 2.0;
    uS = (uSelf + flux[0][k][j-1][i]) / 2.0;
    uN = (uSelf + flux[0][k][j+1][i]) / 2.0;
    uB = (uSelf + flux[0][k-1][j][i]) / 2.0;
    uF = (uSelf + flux[0][k+1][j][i]) / 2.0;

    // prepare v
    vS = (flux[1][k][j-1][i] + flux[1][k][j-1][i+1]) / 2.0;
    vN = (flux[1][k][j][i] + flux[1][k][j][i+1]) / 2.0;

    // prepare w
    wB = (flux[2][k-1][j][i] + flux[2][k-1][j][i+1]) / 2.0;
    wF = (flux[2][k][j][i] + flux[2][k][j][i+1]) / 2.0;

    return 
        (uE * uE - uW * uW) / ctx->mesh->dL[0][0][i] + 
        (vN * uN - vS * uS) / ctx->mesh->dL[0][1][j] + 
        (wF * uF - wB * uB) / ctx->mesh->dL[0][2][k];
}


/** \brief a private kernel for the convection at a v-velocity point in 3D. */
inline PetscReal kernelV(
        NonLinearCtx const * const &ctx, const std::vector<PetscReal***> &flux, 
        const PetscInt &i, const PetscInt &j, const PetscInt &k)
{
    PetscReal   vSelf;
    PetscReal   uW, uE;
    PetscReal   vS, vN, vW, vE, vB, vF;
    PetscReal   wB, wF;

    // prepare self
    vSelf = flux[1][k][j][i];

    // prepare u
    uW = (flux[0][k][j][i-1] + flux[0][k][j+1][i-1]) / 2.0;
    uE = (flux[0][k][j][i] + flux[0][k][j+1][i]) / 2.0;

    // prepare v
    vW = (vSelf + flux[1][k][j][i-1]) / 2.0;
    vE = (vSelf + flux[1][k][j][i+1]) / 2.0;
    vS = (vSelf + flux[1][k][j-1][i]) / 2.0;
    vN = (vSelf + flux[1][k][j+1][i]) / 2.0;
    vB = (vSelf + flux[1][k-1][j][i]) / 2.0;
    vF = (vSelf + flux[1][k+1][j][i]) / 2.0;

    // prepare w
    wB = (flux[2][k-1][j][i] + flux[2][k-1][j+1][i]) / 2.0;
    wF = (flux[2][k][j][i] + flux[2][k][j+1][i]) / 2.0;

    return
        (uE * vE - uW * vW) / ctx->mesh->dL[1][0][i] +
        (vN * vN - vS * vS) / ctx->mesh->dL[1][1][j] +
        (wF * vF - wB * vB) / ctx->mesh->dL[1][2][k];

}


/** \brief a private kernel for the convection at a 3-velocity point in 3D. */
inline PetscReal kernelW(
        NonLinearCtx const * const &ctx, const std::vector<PetscReal***> &flux, 
        const PetscInt &i, const PetscInt &j, const PetscInt &k)
{
    PetscReal   wSelf;
    PetscReal   uW, uE;
    PetscReal   vS, vN;
    PetscReal   wS, wN, wW, wE, wB, wF;

    // prepare self
    wSelf = flux[2][k][j][i];

    // prepare u
    uW = (flux[0][k][j][i-1] + flux[0][k+1][j][i-1]) / 2,0;
    uE = (flux[0][k][j][i] + flux[0][k+1][j][i]) / 2,0;

    // prepare v
    vS = (flux[1][k][j-1][i] + flux[1][k+1][j-1][i]) / 2,0;
    vN = (flux[1][k][j][i] + flux[1][k+1][j][i]) / 2,0;

    // prepare w
    wW = (wSelf + flux[2][k][j][i-1]) / 2.0;
    wE = (wSelf + flux[2][k][j][i+1]) / 2.0;
    wS = (wSelf + flux[2][k][j-1][i]) / 2.0;
    wN = (wSelf + flux[2][k][j+1][i]) / 2.0;
    wB = (wSelf + flux[2][k-1][j][i]) / 2.0;
    wF = (wSelf + flux[2][k+1][j][i]) / 2.0;

    return
        (uE * wE - uW * wW) / ctx->mesh->dL[2][0][i] +
        (vN * wN - vS * wS) / ctx->mesh->dL[2][1][j] +
        (wF * wF - wB * wB) / ctx->mesh->dL[2][2][k];
}
