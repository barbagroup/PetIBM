/***************************************************************************//**
 * \file SimulationParameters.h
 * \author Anush Krishnan (anush@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Definition of the class `SimulationParameters`.
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
# include "FlowDescription.h"
# include "CartesianMesh.h"


class Solutions
{
public:

    PetscInt                dim;

    Vec                     qGlobal;

    Vec                     lambdaGlobal;

    std::string             info;
    
    
    Solutions();

    Solutions(const CartesianMesh &mesh, 
            const types::OutputType &type=types::HDF5);

    ~Solutions();


    PetscErrorCode init(const CartesianMesh &mesh, 
            const types::OutputType &type=types::HDF5);

    PetscErrorCode setOutputFormat(const types::OutputType &type);
    PetscErrorCode setOutputFluxFlag(const PetscBool &flag=PETSC_TRUE);

    PetscErrorCode printInfo() const;

    PetscErrorCode applyIC(const FlowDescription &flow, const Mat &R=nullptr);

    PetscErrorCode write(const std::string &dir, const std::string &name) const;

    PetscErrorCode getVelocity(const Vec &Rinv, Vec &U) const;


protected:

    std::shared_ptr<const MPI_Comm>         comm;

    PetscMPIInt                             mpiSize,
                                            mpiRank;

    std::shared_ptr<const CartesianMesh>    mesh;

    PetscBool                               fluxFlag;

    std::string                             fileExt;
    PetscViewerType                         viewerType;

    PetscErrorCode createInfoString();

private:

};


std::ostream &operator<< (std::ostream &os, const Solutions &soln);
