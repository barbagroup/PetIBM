/**
 * \file misc.cpp
 * \brief Implementations of some miscellaneous functions.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 */

# include <petscsys.h>
# include <petibm/misc.h>
# include <petibm/mesh.h>


namespace petibm
{
namespace misc
{
    PetscErrorCode checkPeriodicBC(
            const type::IntVec2D &bcTypes, type::BoolVec2D &periodic)
    {
        using namespace petibm::type;
        
        PetscFunctionBeginUser;
        
        // dimension
        PetscInt dim = 3;
        
        // initialize periodic bools
        periodic = BoolVec2D(3, BoolVec1D(3, PETSC_FALSE));

        // if a boundary of a field is periodic, the counterpart should be, too
        for(int f=0; f<3; f++)
        {
            for(int d=0; d<3; d++)
            {
                bool p1 = (bcTypes[f][d*2] == int(BCType::PERIODIC));
                bool p2 = (bcTypes[f][d*2+1] == int(BCType::PERIODIC));
                
                if (p1 && (!p2))
                {
                    SETERRQ3(PETSC_COMM_WORLD, PETSC_ERR_ARG_INCOMP,
                            "Boundary %s for velocity field %s is periodic BC, "
                            "but its counterpart, %s, is not!\n",
                            bl2str[BCLoc(d*2)].c_str(), fd2str[Field(f)].c_str(),
                            bl2str[BCLoc(d*2+1)].c_str());
                }
                
                if ((!p1) && p2)
                {
                    SETERRQ3(PETSC_COMM_WORLD, PETSC_ERR_ARG_INCOMP,
                            "Boundary %s for velocity field %s is periodic BC, "
                            "but its counterpart, %s, is not!\n",
                            bl2str[BCLoc(d*2+1)].c_str(), fd2str[Field(f)].c_str(),
                            bl2str[BCLoc(d*2)].c_str());
                }
                
                if (p1 && p2) periodic[f][d] = PETSC_TRUE;
            }
        }
        
        // if all BCs of w velocity are NOBC, it's a 3D simulation
        if (std::all_of(bcTypes[2].begin(), bcTypes[2].end(),
                    [](int t){return t == int(BCType::NOBC);}))
            dim = 2;

        // if a field at a boundary is periodic BC, other fields should be, too
        for(int d=0; d<3; d++)
        {
            PetscInt c = 0;
            for(int f=0; f<3; f++)
            {
                if (periodic[f][d]) c += 1;
            }
            if ((c != 0) && (c != dim))
                SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_ARG_INCOMP,
                        "Not all velocity fields in direction %s are periodic!",
                        type::dir2str[Dir(d)].c_str());
        }

        PetscFunctionReturn(0);
    } // checkPeriodicBC
    

    PetscErrorCode checkBoundaryProc(const DM &da, const type::IntVec1D &n,
            const type::BCLoc &loc, PetscBool &onThisProc)
    {
        PetscFunctionBeginUser;
        
        PetscErrorCode  ierr;

        MPI_Comm    bcMPI;

        switch (loc)
        {
            case type::BCLoc::XMINUS:
                ierr = DMDAGetProcessorSubset(da, DMDA_X, 0, &bcMPI); 
                CHKERRQ(ierr);
                break;
            case type::BCLoc::XPLUS:
                ierr = DMDAGetProcessorSubset(da, DMDA_X, n[0]-1, &bcMPI); 
                CHKERRQ(ierr);
                break;
            case type::BCLoc::YMINUS:
                ierr = DMDAGetProcessorSubset(da, DMDA_Y, 0, &bcMPI); 
                CHKERRQ(ierr);
                break;
            case type::BCLoc::YPLUS:
                ierr = DMDAGetProcessorSubset(da, DMDA_Y, n[1]-1, &bcMPI); 
                CHKERRQ(ierr);
                break;
            case type::BCLoc::ZMINUS:
                ierr = DMDAGetProcessorSubset(da, DMDA_Z, 0, &bcMPI); 
                CHKERRQ(ierr);
                break;
            case type::BCLoc::ZPLUS:
                ierr = DMDAGetProcessorSubset(da, DMDA_Z, n[2]-1, &bcMPI); 
                CHKERRQ(ierr);
                break;
        }

        if (bcMPI != MPI_COMM_NULL)
            onThisProc = PETSC_TRUE;
        else
            onThisProc = PETSC_FALSE;
            PetscFunctionReturn(0);
        } // checkBoundaryProc


    PetscErrorCode getGhostPointList(const type::Mesh &mesh,
            const type::Field &field, const type::BCLoc &loc,
            type::GhostPointsList &points)
    {
        PetscFunctionBeginUser;

        PetscErrorCode      ierr;
        
        // the direction of the normal of this boundary; 1: x, 2: y, 3: z
        PetscInt            axis = int(loc) / 2;

        // the two directions that will be used in the following double loop
        type::IntVec1D      pAxes;
        
        // get pAxes
        ierr = getPerpendAxes(axis, pAxes); CHKERRQ(ierr);

        // alias
        const type::IntVec1D &bg = mesh->bg[field];
        const type::IntVec1D &ed = mesh->ed[field];
        const type::IntVec1D &n = mesh->n[field];
        const type::GhostedVec2D &dL = mesh->dL[field];
        const type::GhostedVec2D &coord = mesh->coord[field];

        for(PetscInt b=bg[pAxes[1]]; b<ed[pAxes[1]]; ++b)
        {
            for(PetscInt a=bg[pAxes[0]]; a<ed[pAxes[0]]; ++a)
            {
                // he stencil of the ghost point
                MatStencil  ghost;
                
                // the stencil of the point associated to the ghost point
                MatStencil  target;

                // targetId is the index of target in packed global Vec
                PetscInt    targetId;

                // ghostId is the index of ghost in a local Vec
                PetscInt    ghostId;
                
                // get the stencils of ghost and target
                ierr = misc::getGhostTargetStencil(
                        n, loc, {a, b}, ghost, target); CHKERRQ(ierr);

                // get packed index of the target
                ierr = mesh->getPackedGlobalIndex(field, target, targetId); CHKERRQ(ierr);
                
                // get the local index of the ghost
                ierr = mesh->getLocalIndex(field, ghost, ghostId); CHKERRQ(ierr);

                // area is the cell surface area used to calculate flux
                // d is the distance between the target and the ghost points
                PetscReal   area, d;

                area = dL[pAxes[0]][a] * dL[pAxes[1]][b];

                if ((int(loc) == 1) || (int(loc) == 3) || (int(loc) == 5))
                    d = coord[axis][n[axis]] - coord[axis][n[axis]-1];
                else if ((int(loc) == 0) || (int(loc) == 2) || (int(loc) == 4))
                    d = coord[axis][0] - coord[axis][-1];

                points[ghost] = {ghostId, target, targetId, area, d, 0.0, 0.0, 0.0};
            }
        }

        PetscFunctionReturn(0);
    } // getGhostPointList
    

    PetscErrorCode getPerpendAxes(const PetscInt &self, type::IntVec1D &pAxes)
    {
        PetscFunctionBeginUser;
        
        pAxes = type::IntVec1D(2, 0);
        
        switch (self)
        {
            case 0:
                pAxes = {1, 2};
                break;
            case 1:
                pAxes = {0, 2};
                break;
            case 2:
                pAxes = {0, 1};
                break;
            default:
                SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                        "Can not recongnize axis %d in the function "
                        "getPerpendAxes!", self);
        }
        
        PetscFunctionReturn(0);
    } // getPerpendAxes
    

    PetscErrorCode getGhostTargetStencil(const type::IntVec1D &n,
            const type::BCLoc &loc, const type::IntVec1D &pIdx,
            MatStencil &ghost, MatStencil &target)
    {
        PetscFunctionBeginUser;
        
        switch (int(loc))
        {
            case 0: // XMINUS
                target = {pIdx[1], pIdx[0], 0, 0};
                ghost = {pIdx[1], pIdx[0], -1, 0};
                break;
            case 1: // XPLUS
                target = {pIdx[1], pIdx[0], n[0]-1, 0};
                ghost = {pIdx[1], pIdx[0], n[0], 0};
                break;
            case 2: // YMINUS
                target = {pIdx[1], 0, pIdx[0], 0};
                ghost = {pIdx[1], -1, pIdx[0], 0};
                break;
            case 3: // YPLUS
                target = {pIdx[1], n[1]-1, pIdx[0], 0};
                ghost = {pIdx[1], n[1], pIdx[0], 0};
                break;
            case 4: // ZMINUS
                target = {0, pIdx[1], pIdx[0], 0};
                ghost = {-1, pIdx[1], pIdx[0], 0};
                break;
            case 5: // ZPLUS
                target = {n[2]-1, pIdx[1], pIdx[0], 0};
                ghost = {n[2], pIdx[1], pIdx[0], 0};
                break;
            default:
                SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                        "Can not recongnize BC location %s in the function "
                        "getGhostTargetStencil!", type::bl2str[loc].c_str());
        }
        
        PetscFunctionReturn(0);
    } // getGhostTargetStencil

} // end of namespace misc
} // end of namespace petibm
