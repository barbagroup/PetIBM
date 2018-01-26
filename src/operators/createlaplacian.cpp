/**
 * \file createlaplacian.cpp
 * \brief Definition of functions creating Laplacian operator.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 */


// STL
# include <map>
# include <numeric>

// here goes PETSc headers
# include <petscmat.h>

// here goes headers from our PetIBM
# include <petibm/mesh.h>
# include <petibm/boundary.h>
# include <petibm/type.h>
# include <petibm/misc.h>


namespace petibm
{
namespace operators
{

/** \brief a private struct for matrix-free constraint mat. */
struct LagrangianCtx
{
    const type::Boundary    bc;
    type::MatrixModifier    modifier;
    
    LagrangianCtx(const type::Boundary &_bc):
        bc(_bc), modifier(bc->dim) {};
};


/** \brief a user-defined MatMult for constraint mat. */
PetscErrorCode LCorrectionMult(Mat mat, Vec x, Vec y);


/** \brief a user-defined MatDestroy for constraint mat. */
PetscErrorCode LCorrectionDestroy(Mat mat);


/** \brief a helper type definition for vector of MatStencil. */
typedef std::vector<MatStencil> StencilVec;


/** \brief a function type for functions returning neighbor stencils. */
typedef std::function<StencilVec(const PetscInt &, 
        const PetscInt &, const PetscInt &)>            GetStencilsFunc;


/**
 * \brief a private function that sets values of a row in Laplacian.
 *
 * \param mesh an instance of CartesianMesh class.
 * \param bc an instance of Boundary class.
 * \param getStencils a function that returns stencils according to dim.
 * \param f the index of field (u=0, v=1, z=2).
 * \param i the x-index of current row in the Cartesian mesh.
 * \param j the y-index of current row in the Cartesian mesh.
 * \param k the z-index of current row in the Cartesian mesh.
 * \param L the Laplacian matrix.
 * \param rowModifiers an object holding information about where in a matrix 
 *                     should be modified.
 *
 * \return PetscErrorCode.
 */
inline PetscErrorCode setRowValues(
        const type::Mesh &mesh,
        const type::Boundary &bc,
        const GetStencilsFunc &getStencils,
        const PetscInt &f, const PetscInt &i,
        const PetscInt &j, const PetscInt &k,
        Mat &L,
        std::map<MatStencil, type::RowModifier> &rowModifiers);


// implementation of petibm::operators::createLaplacian
PetscErrorCode createLaplacian(const type::Mesh &mesh,
                               const type::Boundary &bc, 
                               Mat &L, Mat &LCorrection)
{
    using namespace std::placeholders;

    PetscFunctionBeginUser;

    PetscErrorCode                  ierr;
    
    GetStencilsFunc                 getStencils;

    LagrangianCtx                   *ctx;


    // initialize modifier, regardless what's inside it now
    ctx = new LagrangianCtx(bc);

    // set up getStencils
    if (mesh->dim == 2)
        getStencils = 
            [](const PetscInt &i, const PetscInt &j, const PetscInt &k)
            -> StencilVec {
                return {{k, j, i, 0}, 
                        {k, j, i-1, 0}, {k, j, i+1, 0}, 
                        {k, j-1, i, 0}, {k, j+1, i, 0}};
            };
    else // implies dim == 3
        getStencils = 
            [](const PetscInt &i, const PetscInt &j, const PetscInt &k)
            -> StencilVec {
                return {{k, j, i, 0}, 
                        {k, j, i-1, 0}, {k, j, i+1, 0}, 
                        {k, j-1, i, 0}, {k, j+1, i, 0},
                        {k-1, j, i, 0}, {k+1, j, i, 0}};
            };


    // create matrix
    ierr = DMCreateMatrix(mesh->UPack, &L); CHKERRQ(ierr);
    ierr = MatSetFromOptions(L); CHKERRQ(ierr);
    ierr = MatSetOption(L, MAT_KEEP_NONZERO_PATTERN, PETSC_FALSE); CHKERRQ(ierr);
    ierr = MatSetOption(L, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE); CHKERRQ(ierr);


    // set values, one field each time
    for(PetscInt field=0; field<mesh->dim; ++field)
    {
        // set row values
        // the std::bind originally only binds values, in order to use
        // reference instead of value-copying, we have to use std::ref
        ierr = misc::tripleLoops(
                {mesh->bg[field][2], mesh->ed[field][2]},
                {mesh->bg[field][1], mesh->ed[field][1]},
                {mesh->bg[field][0], mesh->ed[field][0]},
                std::bind(setRowValues, 
                    std::ref(mesh), std::ref(bc), std::ref(getStencils), 
                    std::ref(field), _3, _2, _1, std::ref(L), 
                    std::ref(ctx->modifier[field])));
        CHKERRQ(ierr);
    }


    // temporarily assemble matrix
    ierr = MatAssemblyBegin(L, MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(L, MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);


    // TODO: check if a0 coefficients are already up-to-date.
    for(PetscInt f=0; f<mesh->dim; ++f)
        for(auto &bd: bc->bds[f])
            if (bd->onThisProc)
                for(auto &pt: bd->points)
                {
                    PetscInt    col = pt.second.targetPackedId;
                    PetscReal   value = 
                        ctx->modifier[f][pt.first].coeff * pt.second.a0;

                    ierr = MatSetValue(L, ctx->modifier[f][pt.first].row, 
                            col, value, ADD_VALUES); CHKERRQ(ierr);
                }


    // assemble matrix, an implicit MPI barrier is applied
    ierr = MatAssemblyBegin(L, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(L, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);


    // create a matrix-free constraint matrix for boundary correction
    ierr = MatCreateShell(mesh->comm, mesh->UNLocal, mesh->UNLocal,
            PETSC_DETERMINE, PETSC_DETERMINE, (void *) ctx, &LCorrection); 
    CHKERRQ(ierr);

    ierr = MatShellSetOperation(LCorrection, MATOP_MULT,
            (void(*)(void)) LCorrectionMult); CHKERRQ(ierr);

    ierr = MatShellSetOperation(LCorrection, MATOP_DESTROY,
            (void(*)(void)) LCorrectionDestroy); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // createLaplacian


// implementation of setRowValues
inline PetscErrorCode setRowValues(
        const type::Mesh &mesh,
        const type::Boundary &bc, 
        const GetStencilsFunc &getStencils, 
        const PetscInt &f, const PetscInt &i, 
        const PetscInt &j, const PetscInt &k, Mat &L,
        std::map<MatStencil, type::RowModifier> &rowModifiers)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    StencilVec          stencils = getStencils(i, j, k);

    type::IntVec1D      cols(1+2*mesh->dim, -1);

    type::RealVec1D     values(1+2*mesh->dim, 0.0);


    // get the local index for points in stencils
    for(unsigned int id=0; id<stencils.size(); ++id)
    {
        ierr = mesh->getPackedGlobalIndex(
                f, stencils[id], cols[id]); CHKERRQ(ierr);
    }


    // get the values. One direction each time.
    // Luckily we have dL for ghost points, and the index of -1 is usable in dL,
    // so we don't have to use "if" to find out the boundary points
    for(PetscInt dir=0; dir<mesh->dim; ++dir)
    {
        // determine the index based on the direction
        const PetscInt &self = (dir == 0)? i : (dir == 1)? j : k;

        // alias for cleaner code
        const PetscReal &dLSelf = mesh->dL[f][dir][self];
        const PetscReal &dLNeg = mesh->coord[f][dir][self] - mesh->coord[f][dir][self-1];
        const PetscReal &dLPos = mesh->coord[f][dir][self+1] - mesh->coord[f][dir][self];

        values[dir*2+1] = 1.0 / (dLNeg * dLSelf);
        values[dir*2+2] = 1.0 / (dLPos * dLSelf);
    }

    // update diagonal value
    values[0] = - std::accumulate(values.begin()+1, values.end(), 0.0);

    // set the coefficient for the left point
    ierr = MatSetValues(L, 1, &cols[0], stencils.size(), cols.data(), 
            values.data(), INSERT_VALUES); CHKERRQ(ierr);

    // save the values for boundary ghost points
    for(unsigned int id=0; id<stencils.size(); ++id)
        if (cols[id] == -1) rowModifiers[stencils[id]] = {cols[0], values[id]};

    PetscFunctionReturn(0);
} // setRowValues


// implementation of LCorrectionMult
PetscErrorCode LCorrectionMult(Mat mat, Vec x, Vec y)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    LagrangianCtx       *ctx;

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
    // but there is an implicit MPI barrier in assembly, which we definitely need
    ierr = VecAssemblyBegin(y); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(y); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // LCorrectionMult


// implementation of LCorrectionDestroy
PetscErrorCode LCorrectionDestroy(Mat mat)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;
    LagrangianCtx       *ctx;

    // get the context
    ierr = MatShellGetContext(mat, (void *) &ctx); CHKERRQ(ierr);

    // deallocate the context
    delete ctx;

    PetscFunctionReturn(0);
} // LCorrectionDestroy

} // end of namespace operators
} // end of namespace petibm
