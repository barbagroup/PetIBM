/***************************************************************************//**
 * \file SingleBoundary.cpp
 * \author Anush Krishnan (anus@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Definition of the member functions of class `SingleBoundary`.
 */


// here goes headers from our PetIBM
#include "SingleBoundary.h"


/** \copydoc SingleBoundary::SingleBoundary() */
SingleBoundary::SingleBoundary() = default;


/** \copydoc SingleBoundary::~SingleBoundary(). */
SingleBoundary::~SingleBoundary() = default;


/** \copydoc SingleBoundary::SingleBoundary(const CartesianMesh &, const types::BCLoc &). */
SingleBoundary::SingleBoundary(
        const CartesianMesh &mesh, const types::BCLoc &loc)
{
    init(mesh, loc);
}


/** \copydoc SingleBoundary::init */
PetscErrorCode SingleBoundary::init(
        const CartesianMesh &_mesh, const types::BCLoc &_loc)
{
    using namespace std::placeholders;

    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    // create a shared pointer to mesh; bad practice...
    mesh = std::shared_ptr<const CartesianMesh>(&_mesh, [](const CartesianMesh*){}); 

    // obtain MPI information from CartesianMesh object
    comm = mesh->comm;
    mpiSize = mesh->mpiSize;
    mpiRank = mesh->mpiRank;

    // set dim
    dim = mesh->dim;

    // set the location
    loc = _loc;

    // set onThisProc
    ierr = setProc(); CHKERRQ(ierr);

    // for processes on this boundary, set up IDs and values
    if (onThisProc)
    {
        // initialize the STL vectors
        type.resize(dim); value.resize(dim);
        points.resize(dim);
        updateCoeffsKernel.resize(dim);

        // get the types and values
        for(PetscInt f=0; f<dim; ++f)
        {
            type[f] = (*(mesh->bcInfo))[loc][types::Field(f)].type;
            value[f] = (*(mesh->bcInfo))[loc][types::Field(f)].value;

            ierr = setPoints(f); CHKERRQ(ierr);
        }

        // set corresponding updating functions
        updateCoeffs = std::bind(&SingleBoundary::updateCoeffsTrue, this, _1, _2);
        updateGhosts = std::bind(&SingleBoundary::updateGhostsTrue, this, _1);
    }
    else
    {
        updateCoeffs = 
            [](Solutions &, const PetscReal &)->PetscErrorCode {PetscFunctionReturn(0);};
        updateGhosts = [](Solutions &)->PetscErrorCode {PetscFunctionReturn(0);};
    }

    ierr = MPI_Barrier(*comm); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/** \copydoc SingleBoundary::getProc */
PetscErrorCode SingleBoundary::setProc()
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    MPI_Comm    bcMPI;

    switch (loc)
    {
        case types::BCLoc::XMINUS:
            ierr = DMDAGetProcessorSubset(mesh->da[3], DMDA_X, 0, &bcMPI); 
            CHKERRQ(ierr);
            break;
        case types::BCLoc::XPLUS:
            ierr = DMDAGetProcessorSubset(
                    mesh->da[3], DMDA_X, mesh->n[3][0]-1, &bcMPI); 
            CHKERRQ(ierr);
            break;
        case types::BCLoc::YMINUS:
            ierr = DMDAGetProcessorSubset(mesh->da[3], DMDA_Y, 0, &bcMPI); 
            CHKERRQ(ierr);
            break;
        case types::BCLoc::YPLUS:
            ierr = DMDAGetProcessorSubset(
                    mesh->da[3], DMDA_Y, mesh->n[3][1]-1, &bcMPI); 
            CHKERRQ(ierr);
            break;
        case types::BCLoc::ZMINUS:
            ierr = DMDAGetProcessorSubset(mesh->da[3], DMDA_Z, 0, &bcMPI); 
            CHKERRQ(ierr);
            break;
        case types::BCLoc::ZPLUS:
            ierr = DMDAGetProcessorSubset(
                    mesh->da[3], DMDA_Z, mesh->n[3][2]-1, &bcMPI); 
            CHKERRQ(ierr);
            break;
    }

    if (bcMPI != MPI_COMM_NULL)
        onThisProc = PETSC_TRUE;
    else
        onThisProc = PETSC_FALSE;

    PetscFunctionReturn(0);
}


/** \copydoc SingleBoundary::getIds */
PetscErrorCode SingleBoundary::setPoints(const PetscInt &field)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;
    PetscInt            self, ghost;

    switch (loc)
    {
        case types::BCLoc::XMINUS:
            self = 0; ghost = self - 1;
            normal = -1.0;
            ierr = setPointsX(field, self, ghost); CHKERRQ(ierr);
            break;
        case types::BCLoc::XPLUS:
            self = mesh->n[field][0]-1; ghost = self + 1;
            normal = 1.0;
            ierr = setPointsX(field, self, ghost); CHKERRQ(ierr);
            break;
        case types::BCLoc::YMINUS:
            self = 0; ghost = self - 1;
            normal = -1.0;
            ierr = setPointsY(field, self, ghost); CHKERRQ(ierr);
            break;
        case types::BCLoc::YPLUS:
            self = mesh->n[field][1]-1; ghost = self + 1;
            normal = 1.0;
            ierr = setPointsY(field, self, ghost); CHKERRQ(ierr);
            break;
        case types::BCLoc::ZMINUS:
            self = 0; ghost = self - 1;
            normal = -1.0;
            ierr = setPointsZ(field, self, ghost); CHKERRQ(ierr);
            break;
        case types::BCLoc::ZPLUS:
            self = mesh->n[field][2]-1; ghost = self + 1;
            normal = 1.0;
            ierr = setPointsZ(field, self, ghost); CHKERRQ(ierr);
            break;
    }

    PetscFunctionReturn(0);
}


/** \copydoc SingleBoundary::setPointsX */
PetscErrorCode SingleBoundary::setPointsX(
        const PetscInt &field, const PetscInt &self, const PetscInt &ghost)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    for(PetscInt k=mesh->bg[field][2]; k<mesh->ed[field][2]; ++k)
    {
        for(PetscInt j=mesh->bg[field][1]; j<mesh->ed[field][1]; ++j)
        {
            // selfId is the local id in a Vec of the boundary point
            // ghostId is the local id in a Vec of the ghost point
            // area is the cell surface area used to calculate flux
            // dL is the distance between the boundary and the ghost points
            PetscInt    selfId, ghostId;
            PetscReal   area, dL;

            ierr = DMDAConvertToCell(mesh->da[field], {k, j, self, 0}, &selfId);
            CHKERRQ(ierr);

            ierr = DMDAConvertToCell(mesh->da[field], {k, j, ghost, 0}, &ghostId);
            CHKERRQ(ierr);

            area = mesh->dL[field][1][j] * mesh->dL[field][2][k];

            if (normal == 1)
                // for right boundary, dL is the dx of the last pressure cell
                dL = mesh->dL[3][0][mesh->n[3][0]-1];
            else
                // for left boundary, dL is the dx of the first pressure cell
                dL = mesh->dL[3][0][0];

            points[field].push_back({selfId, ghostId, area, dL, 0.0, 0.0});
        }
    }

    ierr = setFunctions(field, 0); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/** \copydoc SingleBoundary::setPointsY */
PetscErrorCode SingleBoundary::setPointsY(
        const PetscInt &field, const PetscInt &self, const PetscInt &ghost)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    for(PetscInt k=mesh->bg[field][2]; k<mesh->ed[field][2]; ++k)
    {
        for(PetscInt i=mesh->bg[field][0]; i<mesh->ed[field][0]; ++i)
        {
            // selfId is the local id in a Vec of the boundary point
            // ghostId is the local id in a Vec of the ghost point
            // area is the cell surface area used to calculate flux
            // dL is the distance between the boundary and the ghost points
            PetscInt    selfId, ghostId;
            PetscReal   area, dL;

            ierr = DMDAConvertToCell(mesh->da[field], {k, self, i, 0}, &selfId);
            CHKERRQ(ierr);

            ierr = DMDAConvertToCell(mesh->da[field], {k, ghost, i, 0}, &ghostId);
            CHKERRQ(ierr);

            area = mesh->dL[field][0][i] * mesh->dL[field][2][k];

            if (normal == 1)
                // for top boundary, dL is the dy of the last pressure cell
                dL = mesh->dL[3][1][mesh->n[3][1]-1];
            else
                // for bottom boundary, dL is the dx of the first pressure cell
                dL = mesh->dL[3][1][0];

            points[field].push_back({selfId, ghostId, area, dL, 0.0, 0.0});
        }
    }

    ierr = setFunctions(field, 1); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/** \copydoc SingleBoundary::setPointsZ */
PetscErrorCode SingleBoundary::setPointsZ(
        const PetscInt &field, const PetscInt &self, const PetscInt &ghost)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    for(PetscInt j=mesh->bg[field][1]; j<mesh->ed[field][1]; ++j)
    {
        for(PetscInt i=mesh->bg[field][0]; i<mesh->ed[field][0]; ++i)
        {
            // selfId is the local id in a Vec of the boundary point
            // ghostId is the local id in a Vec of the ghost point
            // area is the cell surface area used to calculate flux
            // dL is the distance between the boundary and the ghost points
            PetscInt    selfId, ghostId;
            PetscReal   area, dL;

            ierr = DMDAConvertToCell(mesh->da[field], {self, j, i, 0}, &selfId);
            CHKERRQ(ierr);

            ierr = DMDAConvertToCell(mesh->da[field], {ghost, j, i, 0}, &ghostId);
            CHKERRQ(ierr);

            area = mesh->dL[field][0][i] * mesh->dL[field][1][j];

            if (normal == 1)
                // for front boundary, dL is the dz of the last pressure cell
                dL = mesh->dL[3][2][mesh->n[3][2]-1];
            else
                // for back boundary, dL is the dz of the first pressure cell
                dL = mesh->dL[3][2][0];

            points[field].push_back({selfId, ghostId, area, dL, 0.0, 0.0});
        }
    }

    ierr = setFunctions(field, 2); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/** \copydoc SingleBoundary::setFunctions */
PetscErrorCode SingleBoundary::setFunctions(
        const PetscInt &field, const PetscInt &dir)
{
    using namespace std::placeholders;

    PetscFunctionBeginUser;

    switch (type[field])
    {
        case types::BCType::DIRICHLET:
            if (field == dir)
                updateCoeffsKernel[field] = std::bind(
                        [this] (IdPairs &p, const PetscReal &bc) {
                            p.a1 = bc * p.area;},
                        _1, _4);
            else
                updateCoeffsKernel[field] = std::bind(
                        [this] (IdPairs &p, const PetscReal &bc) {
                            p.a0 = -1.0; 
                            p.a1 = 2.0 * bc * p.area;},
                        _1, _4);
            break;

        case types::BCType::NEUMANN:
            updateCoeffsKernel[field] = std::bind(
                    [this] (IdPairs &p, const PetscReal &bc) {
                        p.a0 = 1.0;
                        p.a1 = this->normal * p.dL * bc * p.area;},
                    _1, _4);
            break;

        case types::BCType::CONVECTIVE:
            if (field == dir)
                updateCoeffsKernel[field] =
                    [this] (IdPairs &p, const PetscReal &bdValue, 
                            const PetscReal &ghValue, const PetscReal &bc,
                            const PetscReal &dt) 
                    { 
                        p.a1 = ghValue - 
                            this->normal * dt * bc * (ghValue - bdValue) / p.dL;
                    };
            else
                updateCoeffsKernel[field] =
                    [this] (IdPairs &p, const PetscReal &bdValue, 
                            const PetscReal &ghValue, const PetscReal &bc,
                            const PetscReal &dt) 
                    {
                        p.a0 = -1.0;
                        p.a1 = ghValue + bdValue - 
                            2.0 * this->normal * dt * bc * (ghValue - bdValue) / p.dL;
                    };
            break;

        case types::BCType::PERIODIC:
            break;
    }

    PetscFunctionReturn(0);
}


/** \copydoc SingleBoundary::updateCoeffsTrue */
PetscErrorCode SingleBoundary::updateCoeffsTrue(
        Solutions &soln, const PetscReal &dt)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    const PetscReal     *arry;

    for(PetscInt f=0; f<dim; ++f)
    {
        ierr = VecGetArrayRead(soln.qLocal[f], &arry); CHKERRQ(ierr);

        for(auto &it: points[f]) 
            updateCoeffsKernel[f](it, arry[it.bcPt], arry[it.ghId], value[f], dt); 

        ierr = VecRestoreArrayRead(soln.qLocal[f], &arry); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}


/** \copydoc SingleBoundary::updateGhostTrue */
PetscErrorCode SingleBoundary::updateGhostsTrue(Solutions &soln)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    PetscReal           *arry;

    for(PetscInt f=0; f<dim; ++f)
    {
        ierr = VecGetArray(soln.qLocal[f], &arry); CHKERRQ(ierr);

        for(auto it: points[f])
            arry[it.ghId] = it.a0 * arry[it.bcPt] + it.a1;

        ierr = VecRestoreArray(soln.qLocal[f], &arry); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}
