/**
 * \file misc.h
 * \brief Prototypes of some miscellaneous functions.
 * \author Anush Krishnan (anus@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \copyright MIT.
 */

# pragma once

// here goes C++ STL
# include <vector>
# include <algorithm>
# include <functional>
# include <cmath>

// here goes PETSc
# include <petscsys.h>
# include <petscdmda.h>

// here goes our own headers
# include <petibm/type.h>
# include <petibm/mesh.h>


/**
 * \defgroup miscModule Miscellaneous functions
 * \brief Collections of the namespaces containing rarely used functions.
 * 
 * This module contains all other public functions that are rarely used by API 
 * users. Four namespaces are categorized into this module:
 * petibm::io, petibm::parser, petibm::delta, and petibm::misc. These functions 
 * are sometimes useful at some specific situations.
 * 
 * Also, user-defined YAML converters are categorized into here.
 * 
 * \ingroup petibm
 */


namespace petibm
{
/**
 * \brief A namespace holding miscellaneous functions.
 * \ingroup miscModule
 * */
namespace misc
{
    /**
     * \brief Check if there is any periodic boundary condition and if these
     *        periodic BCs make sense.
     * \param bcTypes [in] a type::IntVec2D obtained from petibm::parser::parseBCs.
     * \param periodic [out] returned 2D array of `PetscBool` indicating the status of periodic BCs.
     * \return PetscErrorCode
     * 
     * This function checks if the settings of periodic boundary conditions make
     * sense, if there is any. Also, it returns a 2D array of bools that
     * indicates where a specific direction in a specific velocity field is 
     * periodic. For example, assuming `p` is the 2D bools returned by this 
     * function, then if `p[1][1]` is `PETSC_TRUE`, this means the y direction of 
     * v velocity field is periodic.
     * 
     * \ingroup miscModule
     */
    PetscErrorCode checkPeriodicBC(
            const type::IntVec2D &bcTypes, type::BoolVec2D &periodic);
    

    /**
     * \brief Check if a boundary is on this process.
     * \param da [in] the DMDA object where this boundary belongs to.
     * \param n [in] a 1D array of the numbers of the cells in x, y and z directions.
     * \param loc [in] the desired location of this boundary.
     * \param onThisProc [out] a PetscBool indicating if is on this process.
     * \return PetscErrorCode.
     * \ingroup miscModule
     */
    PetscErrorCode checkBoundaryProc(const DM &da, const type::IntVec1D &n,
            const type::BCLoc &loc, PetscBool &onThisProc);


    /**
     * \brief Get a list of ghost points on a desired boundary.
     * \param mesh [in] a type::Mesh object.
     * \param field [in] the desired velocity field.
     * \param loc [in] the location of the desired boundary.
     * \param points [out] return type::GhostPointsList.
     * \return PetscErrorCode.
     * \ingroup miscModule
     */
    PetscErrorCode getGhostPointList(const type::Mesh &mesh,
            const type::Field &field, const type::BCLoc &loc,
            type::GhostPointsList &points);
    
    
    /**
     * \brief An utility to get the perpendicular axes of a desired axis.
     * \param self [in] the desired axis: 0, 1, and 2 represent x, y, and z.
     * \param pAxes [out] returned 1D array of perpendicular axes.
     * \return PetscErrorCode.
     * \ingroup miscModule
     */
    PetscErrorCode getPerpendAxes(const PetscInt &self, type::IntVec1D &pAxes);
    
    
    /**
     * \brief Get the stencils of a desired ghost point and its corresponding 
     *        boundary point.
     * \param n [in] a 1D array of numbers of cells in different directions of desired velocity field.
     * \param loc [in] the location of the desired boundary.
     * \param pIdx [in] a 1D array of the indices of the perpendicular axes.
     * \param ghost [out] the returned stencil of the desired ghost point.
     * \param target [out] the stencil of the corresponding boundary point.
     * \return PetscErrorCode
     * \ingroup miscModule
     * 
     * For example, in a 5x5x5 3D mesh, if we want to get the stencils of the 
     * ghost and corresponding points on the `XPLUS` boundary with `j` and `k` 
     * indices equal to `2` and `3` respectively, then the following arguments 
     * can be used:
     * 
     * ```
     * n = {5, 5, 5};
     * loc = petibm::type::BCLoc::XPLUS;
     * pIdx = {2, 3};
     * ```
     * 
     * Then, the returned ghost and boundary points will be:
     * 
     * ```
     * ghost = {3, 2, 5};
     * target = {3, 2, 4};
     * ```
     * Note, the order of `i`, `j`, and `k` in a `MatStencil` is `{k, j, i}`.
     */
    PetscErrorCode getGhostTargetStencil(const type::IntVec1D &n,
            const type::BCLoc &loc, const type::IntVec1D &pIdx,
            MatStencil &ghost, MatStencil &target);
    
        
    /** 
     * \brief Calculate and return cell sizes of stretched grid in one direction.
     * \param bg [in] start of the stretched region.
     * \param ed [in] end of the stretched region.
     * \param n [in] number of total cells in the stretched region.
     * \param r [in] value of the stretched ratio.
     * \param dL [out] cell sizes.
     * \return PetscErrorCode
     * \ingroup miscModule
     */
    inline PetscErrorCode stretchGrid(
            const PetscReal &bg, const PetscReal &ed, 
            const PetscInt &n, const PetscReal &r, type::RealVec1D &dL)
    {
        PetscFunctionBeginUser;

        dL.resize(n);

        // calculate the size of the first cell
        dL[0] = (ed - bg) * (r - 1.0) / (std::pow(r, n) - 1.0);

        // dL[i] = dL[i-1] * r
        for(auto it=dL.begin()+1; it<dL.end(); ++it) *it = *(it -1) * r;

        PetscFunctionReturn(0);
    }


    /**
     * \brief A helper struct to make looping function easier.
     * \ingroup miscModule
     */
    struct LoopBound {const PetscInt &bg, &ed;};


    /**
     * \brief A helper function to carry out a double loop on a given function.
     * \param bound1 [in] bounds for the outer loop.
     * \param bound2 [in] bounds for the inner loop.
     * \param f [in] the kernel that will be called in the function.
     * \return PetscErrorCode.
     * \ingroup miscModule
     */
    inline PetscErrorCode doubleLoops(
            const LoopBound &bound1, const LoopBound &bound2, 
            const std::function<PetscErrorCode(const PetscInt &, const PetscInt &)> &f)
    {
        PetscFunctionBeginUser;

        PetscErrorCode      ierr;

        for(PetscInt i=bound1.bg; i<bound1.ed; i++)
            for(PetscInt j=bound2.bg; j<bound2.ed; j++)
            {
                ierr = f(i, j); CHKERRQ(ierr);
            }
                

        PetscFunctionReturn(0);
    }


    /**
     * \brief A helper function to carry out a triple loop on a given function.
     * \param bound1 [in] bounds for the most outer loop.
     * \param bound2 [in] bounds for the middle loop.
     * \param bound3 [in] bounds for the most inner loop.
     * \param f [in] the kernel that will be called in the function.
     * \return PetscErrorCode.
     * \ingroup miscModule
     */
    inline PetscErrorCode tripleLoops(
            const LoopBound &bound1, const LoopBound &bound2,
            const LoopBound &bound3,
            const std::function<PetscErrorCode(
                const PetscInt &, const PetscInt &, const PetscInt &)> &f)
    {
        using namespace std::placeholders;

        PetscFunctionBeginUser;

        PetscErrorCode      ierr;

        for(PetscInt i=bound1.bg; i<bound1.ed; ++i)
        {
            ierr = doubleLoops(bound2, bound3, std::bind(f, i, _1, _2));
            CHKERRQ(ierr);
        }

        PetscFunctionReturn(0);
    }

} // end of namespace misc
} // end of namespace petibm
