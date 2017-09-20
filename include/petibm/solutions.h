/***************************************************************************//**
 * \file solutions.h
 * \author Anush Krishnan (anush@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Definition of the class `Solutions`.
 */


# pragma once

// C++ STL
# include <string>
# include <vector>
# include <memory>
# include <functional>

// PETSc
# include <petscsys.h>

// PetIBM
# include <petibm/mesh.h>


namespace petibm
{
namespace solution
{

class Solutions
{
public:

    PetscInt                dim;

    Vec                     UGlobal;

    Vec                     pGlobal;

    std::string             info;
    
    
    Solutions();

    Solutions(const type::Mesh &mesh);

    ~Solutions();


    PetscErrorCode init(const type::Mesh &mesh);

    PetscErrorCode applyIC(const YAML::Node &node);

    PetscErrorCode convert2Velocity(const Mat &Rinv);

    PetscErrorCode convert2Flux(const Mat &R);


protected:

    MPI_Comm                                comm;

    PetscMPIInt                             mpiSize,
                                            mpiRank;

    type::Mesh                              mesh;

    PetscErrorCode createInfoString();

private:

};


} // end of namespace solution
} // end of namespace petibm
