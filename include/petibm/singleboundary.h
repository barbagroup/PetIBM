/***************************************************************************//**
 * \file singleboundary.h
 * \author Anush Krishnan (anus@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Definition of the class `SingleBoundary`.
 */


# pragma once

// STL
# include <string>
# include <memory>
# include <functional>

// here goes PETSc headers
# include <petscsys.h>
# include <petscvec.h>

// here goes headers from our PetIBM
# include "cartesianmesh.h"
# include "solutions.h"
# include "type.h"


namespace petibm
{
namespace boundary
{

class SingleBoundary
{
public:

    PetscInt            dim;

    types::BCLoc        loc;

    PetscBool           onThisProc;

    std::vector<types::BCType>  type;

    types::RealVec1D    value;

    PetscReal           normal;

    std::vector<types::GhostPointsList>   points;



    SingleBoundary();

    SingleBoundary(const CartesianMesh &mesh, const types::BCLoc &loc); 

    ~SingleBoundary();


    PetscErrorCode init(const CartesianMesh &mesh, const types::BCLoc &loc); 

    PetscErrorCode setGhostICs(const Solutions &soln);

    std::function<PetscErrorCode(const Solutions &, const PetscReal &)> updateEqs;

    std::function<PetscErrorCode(const Solutions &)> updateGhostValues;

    std::function<PetscErrorCode(std::vector<Vec> &)> copyValues2LocalVecs;

protected:

    std::shared_ptr<const MPI_Comm>     comm;

    PetscMPIInt                         mpiSize,
                                        mpiRank;


    std::shared_ptr<const CartesianMesh> mesh;


    std::vector<std::function<
        void(types::GhostPointInfo &p, const PetscReal &bdValue, 
                const PetscReal &bc, const PetscReal &dt)>>  updateEqsKernel;


    PetscErrorCode setProc();

    PetscErrorCode setPoints(const PetscInt &field);

    PetscErrorCode setPointsX(
            const PetscInt &field, const PetscInt &self, const PetscInt &ghost);

    PetscErrorCode setPointsY(
            const PetscInt &field, const PetscInt &self, const PetscInt &ghost);

    PetscErrorCode setPointsZ(
            const PetscInt &field, const PetscInt &self, const PetscInt &ghost);

    PetscErrorCode setKernels(
            const PetscInt &field, const PetscInt &dir);

    PetscErrorCode updateEqsTrue(const Solutions &soln, const PetscReal &dt);

    PetscErrorCode updateGhostValuesTrue(const Solutions &soln);

    PetscErrorCode copyValues2LocalVecsTrue(std::vector<Vec> &lclVecs);

private:

};

} // end of namespace boundary
} // end of namespace petibm
