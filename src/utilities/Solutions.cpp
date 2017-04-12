/***************************************************************************//**
 * \file SimulationParameters.h
 * \author Anush Krishnan (anush@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Definition of the class `SimulationParameters`.
 */


// here goes C++ STL
# include <fstream>

// PetIBM
#include "Solutions.h"


using namespace types;

/** \copydoc Solutions::Solutions() */
Solutions::Solutions() = default;


/** \copydoc Solutions::~Solutions() */
Solutions::~Solutions() = default;


/** \copydoc Solutions::Solutions(const CartesianMesh &, const OutputType &) */
Solutions::Solutions(const CartesianMesh &mesh, const OutputType &type)
{
    init(mesh, type);
}


PetscErrorCode Solutions::init(
        const CartesianMesh &inMesh, const OutputType &type)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    // create a shared pointer to mesh; bad practice...
    mesh = std::shared_ptr<const CartesianMesh>(&inMesh, [](const CartesianMesh*){}); 

    // obtain MPI information from CartesianMesh object
    comm = mesh->comm;
    mpiSize = mesh->mpiSize;
    mpiRank = mesh->mpiRank;

    // get a copy of dim
    dim = mesh->dim;

    // create global composite vaector for q
    ierr = DMCreateGlobalVector(mesh->qPack, &qGlobal); CHKERRQ(ierr);

    // create global composite vaector for lambda
    ierr = DMCreateGlobalVector(mesh->lambdaPack, &lambdaGlobal); CHKERRQ(ierr);

    // set output format
    ierr = setOutputFormat(type); CHKERRQ(ierr);

    // set the flag indicating whether to output flux (so far, always true)
    ierr = setOutputFluxFlag(); CHKERRQ(ierr);

    // create info string
    ierr = createInfoString(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


PetscErrorCode Solutions::setOutputFormat(const OutputType &type)
{
    PetscFunctionBeginUser;

    using namespace std::placeholders;

    switch (type)
    {
        case OutputType::Binary:
            fileExt = ".dat";
            viewerType = PETSCVIEWERBINARY;
            break;
        case OutputType::VTK:
            SETERRQ(*comm, 56, "VTK format is not supported anymore!");
        case OutputType::HDF5:
            fileExt = ".h5";
            viewerType = PETSCVIEWERHDF5;
            break;
    }

    PetscFunctionReturn(0);
}


PetscErrorCode Solutions::setOutputFluxFlag(const PetscBool &flag)
{
    PetscFunctionBeginUser;

    fluxFlag = flag;

    PetscFunctionReturn(0);
}


PetscErrorCode Solutions::createInfoString()
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;
    PetscInt            size;

    ierr = VecGetSize(qGlobal, &size); CHKERRQ(ierr);

    std::stringstream   ss;

    ss << std::string(80, '=') << std::endl;
    ss << "Solution Vecters:" << std::endl;
    ss << std::string(80, '=') << std::endl;

    ss << "\tDimension: " << dim << std::endl;
    ss << std::endl;

    ierr = VecGetSize(qGlobal, &size); CHKERRQ(ierr);
    ss << "\tLength of Global Packed q Vector: " << size << std::endl;
    ss << std::endl;

    ierr = VecGetSize(lambdaGlobal, &size); CHKERRQ(ierr);
    ss << "\tLength of Global Packed lambda Vector: " << size << std::endl;
    ss << std::endl;

    info = ss.str();

    PetscFunctionReturn(0);
}


PetscErrorCode Solutions::printInfo() const
{
    PetscFunctionBeginUser;
    
    PetscErrorCode ierr;
    ierr = PetscPrintf(*comm, info.c_str()); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // printInfo


PetscErrorCode Solutions::write(
        const std::string &dir, const std::string &name) const
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    std::string     filePath = dir + "/" + name + fileExt;


    PetscViewer         viewer;

    ierr = PetscViewerCreate(*comm, &viewer); CHKERRQ(ierr);
    ierr = PetscViewerSetType(viewer, viewerType); CHKERRQ(ierr);
    ierr = PetscViewerFileSetMode(viewer, FILE_MODE_WRITE); CHKERRQ(ierr);
    ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);


    // output flux
    if (fluxFlag)
    {
        std::vector<Vec>    qGlobalUnPacked(dim);
        ierr = DMCompositeGetAccessArray(mesh->qPack, qGlobal, dim, 
                nullptr, qGlobalUnPacked.data()); CHKERRQ(ierr);
    
        for(PetscInt i=0; i<dim; ++i)
        {
            std::string     objName = "q" + dir2str[Dir(i)];

            ierr = PetscObjectSetName(
                    PetscObject(qGlobalUnPacked[i]), objName.c_str()); CHKERRQ(ierr);

            ierr = VecView(qGlobalUnPacked[i], viewer); CHKERRQ(ierr);
        }

        ierr = DMCompositeRestoreAccessArray(mesh->qPack, qGlobal, dim,
                nullptr, qGlobalUnPacked.data()); CHKERRQ(ierr);
    }


    // output pressure and surface forces
    PetscInt    numDMs;

    ierr = DMCompositeGetNumberDM(mesh->lambdaPack, &numDMs); CHKERRQ(ierr);

    std::vector<Vec>    lambdaGlobalUnPacked(numDMs);
    ierr = DMCompositeGetAccessArray(mesh->lambdaPack, lambdaGlobal, numDMs, 
            nullptr, lambdaGlobalUnPacked.data()); CHKERRQ(ierr);

    // first we output pressure
    ierr = PetscObjectSetName(
            PetscObject(lambdaGlobalUnPacked[0]), "p"); CHKERRQ(ierr);

    ierr = VecView(lambdaGlobalUnPacked[0], viewer); CHKERRQ(ierr);

    // next, we output surface forces
    for(PetscInt i=1; i<numDMs; ++i)
    {
        std::string     objName = "body" + std::to_string(i);

        ierr = PetscObjectSetName(PetscObject(lambdaGlobalUnPacked[i]), 
                objName.c_str()); CHKERRQ(ierr);

        ierr = VecView(lambdaGlobalUnPacked[i], viewer); CHKERRQ(ierr);
    }

    ierr = DMCompositeRestoreAccessArray(mesh->lambdaPack, lambdaGlobal, numDMs, 
            nullptr, lambdaGlobalUnPacked.data()); CHKERRQ(ierr);

    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // write


PetscErrorCode Solutions::getVelocity(const Vec &RInv, Vec &U) const
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    ierr = VecDuplicate(qGlobal, &U); CHKERRQ(ierr);

    ierr = VecPointwiseMult(U, qGlobal, RInv); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


PetscErrorCode Solutions::applyIC(const FlowDescription &flow, const Mat &R)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;
    std::vector<Vec>    qGlobalUnPacked(dim);

    ierr = DMCompositeGetAccessArray(mesh->qPack, qGlobal, dim,
            nullptr, qGlobalUnPacked.data()); CHKERRQ(ierr);

    for(PetscInt comp=0; comp<dim; ++comp)
    {
        ierr = VecSet(qGlobalUnPacked[comp], flow.IC[comp]); CHKERRQ(ierr);
    }

    ierr = DMCompositeRestoreAccessArray(mesh->qPack, qGlobal, dim,
            nullptr, qGlobalUnPacked.data()); CHKERRQ(ierr);

    if (R)
    {
        Vec     y;
        ierr = VecDuplicate(qGlobal, &y); CHKERRQ(ierr);
        ierr = MatMult(R, qGlobal, y); CHKERRQ(ierr);
        ierr = VecSwap(qGlobal, y); CHKERRQ(ierr);
        ierr = VecDestroy(&y); CHKERRQ(ierr);
    }

    ierr = VecSet(lambdaGlobal, 0.0); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


std::ostream &operator<< (std::ostream &os, const Solutions &soln)
{
    os << soln.info;
    return os;
}
