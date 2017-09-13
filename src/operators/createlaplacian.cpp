/***************************************************************************//**
 * \file createlaplacian.cpp
 * \author Anush Krishnan (anus@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Definition of functions regarding to creating Laplacian operators.
 */


// STL
# include <map>
# include <numeric>

// here goes PETSc headers
# include <petscmat.h>

// here goes headers from our PetIBM
# include "petibm/cartesianmesh.h"
# include "petibm/boundary.h"
# include "petibm/type.h"
# include "petibm/misc.h"


namespace petibm
{
namespace operators
{

/** \brief a private struct for matrix-free constraint mat. */
struct LagrangianCtx
{
    std::shared_ptr<const utilities::Boundary>     bc;
    utilities::types::MatrixModifier               modifier;
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
 * \brief a private function that set values of a row in Laplacian.
 *
 * \param mesh an instance of CartesianMesh class.
 * \param bc an instance of Boundary class.
 * \param getStencils a function that returns stencils according to dim.
 * \param f the index of field (u=0, v=1, z=2).
 * \param i the x-index of current row in the cartesian mesh.
 * \param j the y-index of current row in the cartesian mesh.
 * \param k the z-index of current row in the cartesian mesh.
 * \param L the Laplacian matrix.
 *
 * \return PetscErrorCode.
 */
inline PetscErrorCode setRowValues(
		const utilities::CartesianMesh &mesh,
		const utilities::Boundary &bc,
		const GetStencilsFunc &getStencils,
		const PetscInt &f, const PetscInt &i,
		const PetscInt &j, const PetscInt &k,
		Mat &L,
		std::map<MatStencil, utilities::types::RowModifier> &rowModifiers);


/** \copydoc createLaplacian. */
PetscErrorCode createLaplacian(const utilities::CartesianMesh &mesh,
                               const utilities::Boundary &bc, 
                               Mat &L, Mat &LCorrection)
{
    using namespace std::placeholders;

    PetscFunctionBeginUser;

    PetscErrorCode                  ierr;
    
    GetStencilsFunc                 getStencils;

    LagrangianCtx                   *ctx;


    // initialize modifier, regardless what's inside it now
    ctx = new LagrangianCtx;
    ctx->modifier = utilities::types::MatrixModifier(mesh.dim);
    ctx->bc = std::shared_ptr<const utilities::Boundary>(
    		&bc, [](const utilities::Boundary*){});


    // set up getStencils
    if (mesh.dim == 2)
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
    ierr = DMCreateMatrix(mesh.UPack, &L); CHKERRQ(ierr);
    ierr = MatSetFromOptions(L); CHKERRQ(ierr);
    ierr = MatSetOption(L, MAT_KEEP_NONZERO_PATTERN, PETSC_FALSE); CHKERRQ(ierr);
    ierr = MatSetOption(L, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE); CHKERRQ(ierr);


    // set values, one field each time
    for(PetscInt field=0; field<mesh.dim; ++field)
    {
        // set row values
        // the std::bind originally only binds values, in order to use
        // reference instead of value-copying, we have to use std::ref
        ierr = utilities::misc::tripleLoops(
                {mesh.bg[field][2], mesh.ed[field][2]},
                {mesh.bg[field][1], mesh.ed[field][1]},
                {mesh.bg[field][0], mesh.ed[field][0]},
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
    for(auto &bd: bc.bds)
        if (bd.onThisProc)
            for(PetscInt f=0; f<mesh.dim; ++f)
                for(auto &pt: bd.points[f])
                {
                    PetscInt    col = pt.second.targetPackedId;
                    PetscReal   value = 
                        ctx->modifier[f][pt.first].coeff * pt.second.a0;

                    ierr = MatSetValue(L, ctx->modifier[f][pt.first].row, 
                            col, value, ADD_VALUES); CHKERRQ(ierr);
                }


    // assemble matrix, an implicit mpi barrier is applied
    ierr = MatAssemblyBegin(L, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(L, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);


    // create a matrix-free constraint matrix for boundary correction
    ierr = MatCreateShell(*mesh.comm, mesh.UNLocal, mesh.UNLocal,
            PETSC_DETERMINE, PETSC_DETERMINE, (void *) ctx, &LCorrection); 
    CHKERRQ(ierr);

    ierr = MatShellSetOperation(LCorrection, MATOP_MULT,
            (void(*)(void)) LCorrectionMult); CHKERRQ(ierr);

    ierr = MatShellSetOperation(LCorrection, MATOP_DESTROY,
            (void(*)(void)) LCorrectionDestroy); CHKERRQ(ierr);


    // see if users want to view this matrix through cmd
    ierr = PetscObjectViewFromOptions(
            (PetscObject) L, nullptr, "-L_mat_view"); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/** \copydoc setRowValues. */
inline PetscErrorCode setRowValues(
		const utilities::CartesianMesh &mesh,
		const utilities::Boundary &bc, 
		const GetStencilsFunc &getStencils, 
		const PetscInt &f, const PetscInt &i, 
		const PetscInt &j, const PetscInt &k, Mat &L,
		std::map<MatStencil, utilities::types::RowModifier> &rowModifiers)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    StencilVec          stencils = getStencils(i, j, k);

    utilities::types::IntVec1D     lclIds(1+2*mesh.dim, -1),
                                   cols(1+2*mesh.dim, -1);

    utilities::types::RealVec1D    values(1+2*mesh.dim, 0.0);


    // get the local index for points in stencils
    for(unsigned int id=0; id<stencils.size(); ++id)
    {
        ierr = DMDAConvertToCell(
                mesh.da[f], stencils[id], &lclIds[id]); CHKERRQ(ierr);
    }

    // get the global index
    ierr = ISLocalToGlobalMappingApply(
            mesh.UMapping[f], stencils.size(), lclIds.data(), cols.data()); 
    CHKERRQ(ierr);


    // get the values. One direction each time.
    // Luckly we have dL for ghost points, and the index of -1 is usable in dL,
    // so we don't have to use "if" to find out the boundary points
    for(PetscInt dir=0; dir<mesh.dim; ++dir)
    {
        // determine the index based on the direction
        const PetscInt &self = (dir == 0)? i : (dir == 1)? j : k;

        // alias for cleaner code
        const PetscReal &dLSelf = mesh.dL[f][dir][self];
        const PetscReal &dLNeg = mesh.coord[f][dir][self] - mesh.coord[f][dir][self-1];
        const PetscReal &dLPos = mesh.coord[f][dir][self+1] - mesh.coord[f][dir][self];

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
}


/** \copydoc LCorrectionMult. */
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
    for(auto &bd: ctx->bc->bds)
        if (bd.onThisProc)
            for(PetscInt f=0; f<ctx->bc->dim; ++f)
                for(auto &pt: bd.points[f])
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


/** \copydoc LCorrectionDestroy. */
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
}

} // end of namespace operators
} // end of namespace petibm
