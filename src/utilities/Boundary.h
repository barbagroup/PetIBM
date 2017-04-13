/***************************************************************************//**
 * \file Boundary.h
 * \author Anush Krishnan (anus@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Definition of the class `Boundary`.
 */


# pragma once

// STL
# include <vector>

// here goes headers from our PetIBM
# include "CartesianMesh.h"
# include "SingleBoundary.h"
# include "Solutions.h"
# include "types.h"


class Boundary
{
public:

    PetscInt                        dim;

    std::vector<SingleBoundary>     bds;


    Boundary();

    Boundary(const CartesianMesh &mesh);

    ~Boundary();


    PetscErrorCode init(const CartesianMesh &mesh);

    PetscErrorCode setGhostICs(const Solutions &soln);

    PetscErrorCode updateEqs(const Solutions &soln, const PetscReal &dt);

    PetscErrorCode updateGhostValues(const Solutions &soln);

    PetscErrorCode copyValues2LocalVecs(std::vector<Vec> &lclVecs) const;

protected:

    std::shared_ptr<const MPI_Comm>     comm;

    PetscMPIInt                         mpiSize,
                                        mpiRank;


    std::shared_ptr<const CartesianMesh> mesh;

private:

};
