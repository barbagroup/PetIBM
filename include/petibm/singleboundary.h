/***************************************************************************//**
 * \file singleboundary.h
 * \author Anush Krishnan (anus@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Definition of the class `SingleBoundary`.
 */


# pragma once

// C++ STL
# include <memory>

// here goes PETSc headers
# include <petscsys.h>
# include <petscvec.h>

// here goes headers from our PetIBM
# include <petibm/type.h>
# include <petibm/mesh.h>


namespace petibm
{
namespace boundary
{

class SingleBoundaryBase
{
public:

    PetscInt            dim;

    type::BCLoc         loc;
    
    type::Field         field;

    PetscReal           value;

    PetscReal           normal;

    type::GhostPointsList   points;

    PetscBool           onThisProc;



    SingleBoundaryBase() = default;

    SingleBoundaryBase(const type::Mesh &mesh, const type::BCLoc &loc, 
            const type::Field &field, const PetscReal &value); 

    virtual ~SingleBoundaryBase() = default;


    PetscErrorCode setGhostICs(const Vec &vec);

    PetscErrorCode updateEqs(const Vec &vec, const PetscReal &dt);

    PetscErrorCode updateGhostValues(const Vec &vec);

    PetscErrorCode copyValues2LocalVec(Vec &lclVec);

protected:

    PetscErrorCode init(const type::Mesh &mesh, const type::BCLoc &loc, 
            const type::Field &field, const PetscReal &bcValue); 

    virtual PetscErrorCode setGhostICsKernel(
            const PetscReal &targetValue, type::GhostPointInfo &p) = 0;

    virtual PetscErrorCode updateEqsKernel(const PetscReal &targetValue,
            const PetscReal &dt, type::GhostPointInfo &p) = 0;


    MPI_Comm        comm;

    PetscMPIInt     mpiSize;
    
    PetscMPIInt     mpiRank;


    type::Mesh      mesh;

};

} // end of namespace boundary


namespace type
{
    typedef std::shared_ptr<boundary::SingleBoundaryBase> SingleBoundary;
}

namespace boundary
{
    PetscErrorCode createSingleBoundary(
            const type::Mesh &mesh, const type::BCLoc &loc, 
            const type::Field &field, const PetscReal &value,
            type::SingleBoundary &singleBd);
}

} // end of namespace petibm
