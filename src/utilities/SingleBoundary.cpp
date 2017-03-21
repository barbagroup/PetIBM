/***************************************************************************//**
 * \file SingleBoundary.cpp
 * \author Anush Krishnan (anus@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Definition of the member functions of class `SingleBoundary`.
 */


// here goes headers from our PetIBM
#include "SingleBoundary.h"


using namespace types;


/** \copydoc SingleBoundary::SingleBoundary() */
SingleBoundary::SingleBoundary() {}

/** \copydoc SingleBoundary::~SingleBoundary(). */
SingleBoundary::~SingleBoundary() {}

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
        updateCoeffsFuncs.resize(dim);

        // get the types and values
        for(PetscInt f=0; f<dim; ++f)
        {
            type[f] = (*(mesh->bcInfo))[loc][Field(f)].type;
            value[f] = (*(mesh->bcInfo))[loc][Field(f)].value;

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
        case BCLoc::XMINUS:
            ierr = DMDAGetProcessorSubset(mesh->da[3], DMDA_X, 0, &bcMPI); 
            CHKERRQ(ierr);
            break;
        case BCLoc::XPLUS:
            ierr = DMDAGetProcessorSubset(
                    mesh->da[3], DMDA_X, mesh->n[3][0]-1, &bcMPI); 
            CHKERRQ(ierr);
            break;
        case BCLoc::YMINUS:
            ierr = DMDAGetProcessorSubset(mesh->da[3], DMDA_Y, 0, &bcMPI); 
            CHKERRQ(ierr);
            break;
        case BCLoc::YPLUS:
            ierr = DMDAGetProcessorSubset(
                    mesh->da[3], DMDA_Y, mesh->n[3][1]-1, &bcMPI); 
            CHKERRQ(ierr);
            break;
        case BCLoc::ZMINUS:
            ierr = DMDAGetProcessorSubset(mesh->da[3], DMDA_Z, 0, &bcMPI); 
            CHKERRQ(ierr);
            break;
        case BCLoc::ZPLUS:
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
        case BCLoc::XMINUS:
            self = 0; ghost = self - 1;
            normal = -1.0;
            ierr = setPointsX(field, self, ghost); CHKERRQ(ierr);
            break;
        case BCLoc::XPLUS:
            self = mesh->n[field][0]-1; ghost = self + 1;
            normal = 1.0;
            ierr = setPointsX(field, self, ghost); CHKERRQ(ierr);
            break;
        case BCLoc::YMINUS:
            self = 0; ghost = self - 1;
            normal = -1.0;
            ierr = setPointsY(field, self, ghost); CHKERRQ(ierr);
            break;
        case BCLoc::YPLUS:
            self = mesh->n[field][1]-1; ghost = self + 1;
            normal = 1.0;
            ierr = setPointsY(field, self, ghost); CHKERRQ(ierr);
            break;
        case BCLoc::ZMINUS:
            self = 0; ghost = self - 1;
            normal = -1.0;
            ierr = setPointsZ(field, self, ghost); CHKERRQ(ierr);
            break;
        case BCLoc::ZPLUS:
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

            points.push_back({selfId, ghostId, area, dL, 0.0, 0.0});
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

            points.push_back({selfId, ghostId, area, dL, 0.0, 0.0});
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

            points.push_back({selfId, ghostId, area, dL, 0.0, 0.0});
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
        case BCType::DIRICHLET:
            if (field == dir)
                updateCoeffsFuncs[field] = std::bind(
                        &SingleBoundary::updateCoeffsDirichletSameDir, 
                        this, _1, _2, _3);
            else
                updateCoeffsFuncs[field] = std::bind(
                        &SingleBoundary::updateCoeffsDirichletDiffDir, 
                        this, _1, _2, _3);
            break;
        case BCType::NEUMANN:
            if (field == dir)
                updateCoeffsFuncs[field] = std::bind(
                        &SingleBoundary::updateCoeffsNeumannSameDir, 
                        this, _1, _2, _3);
            else
                updateCoeffsFuncs[field] = std::bind(
                        &SingleBoundary::updateCoeffsNeumannDiffDir, 
                        this, _1, _2, _3);
            break;
        case BCType::CONVECTIVE:
            if (field == dir)
                updateCoeffsFuncs[field] = std::bind(
                        &SingleBoundary::updateCoeffsConvectiveSameDir, 
                        this, _1, _2, _3);
            else
                updateCoeffsFuncs[field] = std::bind(
                        &SingleBoundary::updateCoeffsConvectiveDiffDir, 
                        this, _1, _2, _3);
            break;
        case BCType::PERIODIC:
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

        ierr = updateCoeffsFuncs[f](value[f], arry, dt); CHKERRQ(ierr);

        ierr = VecRestoreArrayRead(soln.qLocal[f], &arry); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}


/** \copydoc SingleBoundary::updateCoeffsDirichletSameDir */
PetscErrorCode SingleBoundary::updateCoeffsDirichletSameDir(
        const PetscReal &v, const PetscReal *&arry, const PetscReal &dt)
{
    PetscFunctionBeginUser;
    for(auto &it: points) it.a1 = v * it.area;
    PetscFunctionReturn(0);
}


/** \copydoc SingleBoundary::updateCoeffsDirichletDiffDir */
PetscErrorCode SingleBoundary::updateCoeffsDirichletDiffDir(
        const PetscReal &v, const PetscReal *&arry, const PetscReal &dt)
{
    PetscFunctionBeginUser;
    for(auto &it: points)
    {
        it.a0 = -1.0;
        it.a1 = 2.0 * v * it.area;
    }
    PetscFunctionReturn(0);
}


/** \copydoc SingleBoundary::updateCoeffsNeumannSameDir */
PetscErrorCode SingleBoundary::updateCoeffsNeumannSameDir(
        const PetscReal &v, const PetscReal *&arry, const PetscReal &dt)
{
    PetscFunctionBeginUser;
    for(auto &it: points)
    {
        it.a0 = 1.0;
        it.a1 = normal * it.dL * v * it.area;
    }
    PetscFunctionReturn(0);
}


/** \copydoc SingleBoundary::updateCoeffsNeumannDiffDir */
PetscErrorCode SingleBoundary::updateCoeffsNeumannDiffDir(
        const PetscReal &v, const PetscReal *&arry, const PetscReal &dt)
{
    PetscFunctionBeginUser;
    for(auto &it: points)
    {
        it.a0 = 1.0;
        it.a1 = normal * it.dL * v * it.area;
    }
    PetscFunctionReturn(0);
}


/** \copydoc SingleBoundary::updateCoeffsConvectiveSameDir */
PetscErrorCode SingleBoundary::updateCoeffsConvectiveSameDir(
        const PetscReal &v, const PetscReal *&arry, const PetscReal &dt)
{
    PetscFunctionBeginUser;
    for(auto &it: points)
    {
        it.a1 = arry[it.ghId] - 
            normal * dt * v * (arry[it.ghId] - arry[it.bcPt]) / it.dL;
    }
    PetscFunctionReturn(0);
}


/** \copydoc SingleBoundary::updateCoeffsConvectiveDiffDir */
PetscErrorCode SingleBoundary::updateCoeffsConvectiveDiffDir(
        const PetscReal &v, const PetscReal *&arry, const PetscReal &dt)
{
    PetscFunctionBeginUser;
    for(auto &it: points)
    {
        it.a0 = -1.0;
        it.a1 = arry[it.ghId] + arry[it.bcPt] - 
            2.0 * normal * dt * v * (arry[it.ghId] - arry[it.bcPt]) / it.dL;
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

        for(auto it: points) arry[it.ghId] = it.a0 * arry[it.bcPt] + it.a1;

        ierr = VecRestoreArray(soln.qLocal[f], &arry); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}
