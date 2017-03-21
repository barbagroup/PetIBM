/***************************************************************************//**
 * \file SingleBoundary.h
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
# include "CartesianMesh.h"
# include "Solutions.h"
# include "types.h"


class SingleBoundary
{
public:

    PetscInt            dim;

    types::BCLoc        loc;

    PetscBool           onThisProc;

    std::vector<types::BCType>  type;

    types::RealVec1D    value;

    PetscReal           normal;



    SingleBoundary();
    SingleBoundary(const CartesianMesh &mesh, const types::BCLoc &loc); 

    ~SingleBoundary();


    PetscErrorCode init(const CartesianMesh &mesh, const types::BCLoc &loc); 

    std::function<PetscErrorCode(Solutions &, const PetscReal &)> updateCoeffs;

    std::function<PetscErrorCode(Solutions &)> updateGhosts;

protected:

    std::shared_ptr<const MPI_Comm>     comm;

    PetscMPIInt                         mpiSize,
                                        mpiRank;


    std::shared_ptr<const CartesianMesh> mesh;


    struct IdPairs {PetscInt bcPt; PetscInt ghId; 
        PetscReal area; PetscReal dL; PetscReal a0; PetscReal a1;};

    std::vector<IdPairs>    points;



    PetscErrorCode setProc();

    PetscErrorCode setPoints(const PetscInt &field);

    PetscErrorCode setPointsX(
            const PetscInt &field, const PetscInt &self, const PetscInt &ghost);
    PetscErrorCode setPointsY(
            const PetscInt &field, const PetscInt &self, const PetscInt &ghost);
    PetscErrorCode setPointsZ(
            const PetscInt &field, const PetscInt &self, const PetscInt &ghost);

    PetscErrorCode setFunctions(
            const PetscInt &field, const PetscInt &dir);


    PetscErrorCode updateCoeffsTrue(Solutions &soln, const PetscReal &dt);

    std::vector<std::function<PetscErrorCode(const PetscReal &v, 
            const PetscReal *&, const PetscReal &)>>        updateCoeffsFuncs;

    PetscErrorCode updateCoeffsDirichletSameDir(
            const PetscReal &v, const PetscReal * &arry, const PetscReal &dt);
    PetscErrorCode updateCoeffsDirichletDiffDir(
            const PetscReal &v, const PetscReal * &arry, const PetscReal &dt);
    PetscErrorCode updateCoeffsNeumannSameDir(
            const PetscReal &v, const PetscReal * &arry, const PetscReal &dt);
    PetscErrorCode updateCoeffsNeumannDiffDir(
            const PetscReal &v, const PetscReal * &arry, const PetscReal &dt);
    PetscErrorCode updateCoeffsConvectiveSameDir(
            const PetscReal &v, const PetscReal * &arry, const PetscReal &dt);
    PetscErrorCode updateCoeffsConvectiveDiffDir(
            const PetscReal &v, const PetscReal * &arry, const PetscReal &dt);



    PetscErrorCode updateGhostsTrue(Solutions &soln);

private:

};
