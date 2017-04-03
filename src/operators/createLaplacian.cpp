/***************************************************************************//**
 * \file createLaplacian.cpp
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
# include "CartesianMesh.h"
# include "Boundary.h"
# include "types.h"
# include "misc.h"


/** \brief a helper type definition for vector of MatStencil. */
typedef std::vector<MatStencil> StencilVec;


/** \brief a function type for functions returning neighbor stencils. */
typedef std::function<StencilVec(const PetscInt &, 
        const PetscInt &, const PetscInt &)>            GetStencilsFunc;


typedef std::vector<std::map<PetscInt, PetscReal>>       GhostCoeffType;


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
        const CartesianMesh &mesh, const Boundary &bc, 
        const GetStencilsFunc &getStencils, 
        const PetscInt &f, const PetscInt &i, 
        const PetscInt &j, const PetscInt &k, Mat &L,
        std::map<PetscInt, PetscReal> &bcCoeffs);


/** \copydoc createLaplacian. */
PetscErrorCode createLaplacian(
        const CartesianMesh &mesh, const Boundary &bc, 
        Mat &L, GhostCoeffType &ghCoeff)
{
    using namespace std::placeholders;

    PetscFunctionBeginUser;

    PetscErrorCode                  ierr;
    
    GetStencilsFunc                 getStencils;


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
    ierr = DMCreateMatrix(mesh.qPack, &L); CHKERRQ(ierr);
    ierr = MatSetFromOptions(L); CHKERRQ(ierr);
    ierr = MatSetOption(L, MAT_KEEP_NONZERO_PATTERN, PETSC_FALSE); CHKERRQ(ierr);
    ierr = MatSetOption(L, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE); CHKERRQ(ierr);


    // set values, one field each time
    for(PetscInt field=0; field<mesh.dim; ++field)
    {
        // set row values
        // the std::bind originally only binds values, in order to use
        // reference instead of value-copying, we have to use std::ref
        ierr = misc::tripleLoops(
                {mesh.bg[field][2], mesh.ed[field][2]},
                {mesh.bg[field][1], mesh.ed[field][1]},
                {mesh.bg[field][0], mesh.ed[field][0]},
                std::bind(setRowValues, 
                    std::ref(mesh), std::ref(bc), std::ref(getStencils), 
                    std::ref(field), _3, _2, _1, std::ref(L), 
                    std::ref(ghCoeff[field])));
        CHKERRQ(ierr);
    }


    // TODO: change diagonal entries for BCs.
    for(PetscInt f=0; f<mesh.dim; ++f)
    {
    }


    // assemble matrix
    ierr = MatAssemblyBegin(L, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(L, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);


    // see if users want to view this matrix through cmd
    ierr = PetscObjectViewFromOptions(
            (PetscObject) L, nullptr, "-L_mat_view"); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/** \copydoc setRowValues. */
inline PetscErrorCode setRowValues(
        const CartesianMesh &mesh, const Boundary &bc, 
        const GetStencilsFunc &getStencils, 
        const PetscInt &f, const PetscInt &i, 
        const PetscInt &j, const PetscInt &k, Mat &L,
        std::map<PetscInt, PetscReal> &bcCoeffs)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    StencilVec          stencils = getStencils(i, j, k);

    types::IntVec1D     lclIds(1+2*mesh.dim, -1),
                        cols(1+2*mesh.dim, -1);

    types::RealVec1D    values(1+2*mesh.dim, 0.0);


    // get the local index for points in stencils
    for(PetscInt id=0; id<stencils.size(); ++id)
    {
        ierr = DMDAConvertToCell(
                mesh.da[f], stencils[id], &lclIds[id]); CHKERRQ(ierr);
    }

    // get the global index
    ierr = ISLocalToGlobalMappingApply(
            mesh.qMapping[f], stencils.size(), lclIds.data(), cols.data()); 
    CHKERRQ(ierr);


    // get the values. One direction each time.
    for(PetscInt dir=0; dir<mesh.dim; ++dir)
    {
        const types::BCType &type = 
            (*mesh.bcInfo)[types::BCLoc(dir*2)][types::Field(f)].type;

        const types::RealVec1D &dL = mesh.dL[f][dir];

        const PetscInt &self = (dir == 0)? i : (dir == 1)? j : k;

        const PetscInt &n = dL.size();

        PetscReal   dLNeg = (self > 0) ? dL[self-1] :
            (type == types::PERIODIC) ? dL.back() : mesh.dL[3][dir].front();

        PetscReal   dLPos = (self < (n - 1)) ? dL[self+1] :
            (type == types::PERIODIC) ? dL.front() : mesh.dL[3][dir].back();

        values[dir*2+1] = 1.0 / (0.5 * (dLNeg + dL[self]) * dL[self]);
        values[dir*2+2] = 1.0 / (0.5 * (dL[self] + dLPos) * dL[self]);
    }

    // update diagonal value
    values[0] = - std::accumulate(values.begin()+1, values.end(), 0.0);

    // set the coefficient for the left point
    ierr = MatSetValues(L, 1, &cols[0], stencils.size(), cols.data(), 
            values.data(), INSERT_VALUES); CHKERRQ(ierr);

    // save the values for boundary ghost points
    for(PetscInt id=0; id<stencils.size(); ++id)
        if (cols[id] == -1) bcCoeffs[lclIds[id]] = values[id];

    PetscFunctionReturn(0);
}
