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
    ierr = DMCreateGlobalVector(mesh->UPack, &UGlobal); CHKERRQ(ierr);

    // create global composite vaector for pressure
    ierr = DMCreateGlobalVector(mesh->da[3], &pGlobal); CHKERRQ(ierr);

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

    ierr = VecGetSize(UGlobal, &size); CHKERRQ(ierr);

    std::stringstream   ss;

    ss << std::string(80, '=') << std::endl;
    ss << "Solution Vecters:" << std::endl;
    ss << std::string(80, '=') << std::endl;

    ss << "\tDimension: " << dim << std::endl;
    ss << std::endl;

    ierr = VecGetSize(UGlobal, &size); CHKERRQ(ierr);
    ss << "\tLength of Global Packed Velocity Vector: " << size << std::endl;
    ss << std::endl;

    ierr = VecGetSize(pGlobal, &size); CHKERRQ(ierr);
    ss << "\tLength of Global Pressure Vector: " << size << std::endl;
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

    std::string         filePath = dir + "/" + name + fileExt;

    std::vector<Vec>    UGlobalUnPacked(dim);

    PetscViewer         viewer;


    // create viewer
    ierr = PetscViewerCreate(*comm, &viewer); CHKERRQ(ierr);
    ierr = PetscViewerSetType(viewer, viewerType); CHKERRQ(ierr);
    ierr = PetscViewerFileSetMode(viewer, FILE_MODE_WRITE); CHKERRQ(ierr);
    ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);


    // output velocity
    ierr = DMCompositeGetAccessArray(mesh->UPack, UGlobal, dim, 
            nullptr, UGlobalUnPacked.data()); CHKERRQ(ierr);

    for(PetscInt i=0; i<dim; ++i)
    {
        std::string     objName = types::fd2str[types::Field(i)];

        ierr = PetscObjectSetName(
                PetscObject(UGlobalUnPacked[i]), objName.c_str()); CHKERRQ(ierr);

        ierr = VecView(UGlobalUnPacked[i], viewer); CHKERRQ(ierr);
    }

    ierr = DMCompositeRestoreAccessArray(mesh->UPack, UGlobal, dim,
            nullptr, UGlobalUnPacked.data()); CHKERRQ(ierr);

    // output flux
    if (fluxFlag)
    {
        // TODO: thinks about if outputing flux is necessary
        ierr = PetscPrintf(*comm, "The outputFlux has been set to true, "
                "while current version of PetIBM does not support it yet. "
                "No flux will be output.\n"); CHKERRQ(ierr);
    }


    // output pressure
    ierr = PetscObjectSetName(
            PetscObject(pGlobal), "p"); CHKERRQ(ierr);

    ierr = VecView(pGlobal, viewer); CHKERRQ(ierr);


    // destroy PetscViewer
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // write


PetscErrorCode Solutions::read(
        const std::string &dir, const std::string &name) const
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    std::string         filePath = dir + "/" + name + fileExt;

    std::vector<Vec>    UGlobalUnPacked(dim);

    PetscViewer         viewer;


    // create viewer
    ierr = PetscViewerCreate(*comm, &viewer); CHKERRQ(ierr);
    ierr = PetscViewerSetType(viewer, viewerType); CHKERRQ(ierr);
    ierr = PetscViewerFileSetMode(viewer, FILE_MODE_READ); CHKERRQ(ierr);
    ierr = PetscViewerFileSetName(viewer, filePath.c_str()); CHKERRQ(ierr);


    // read velocity
    ierr = DMCompositeGetAccessArray(mesh->UPack, UGlobal, dim, 
            nullptr, UGlobalUnPacked.data()); CHKERRQ(ierr);

    for(PetscInt i=0; i<dim; ++i)
    {
        std::string     objName = types::fd2str[types::Field(i)];

        ierr = PetscObjectSetName(
                PetscObject(UGlobalUnPacked[i]), objName.c_str()); CHKERRQ(ierr);

        ierr = VecLoad(UGlobalUnPacked[i], viewer); CHKERRQ(ierr);
    }

    ierr = DMCompositeRestoreAccessArray(mesh->UPack, UGlobal, dim,
            nullptr, UGlobalUnPacked.data()); CHKERRQ(ierr);


    // read pressure
    ierr = PetscObjectSetName(
            PetscObject(pGlobal), "p"); CHKERRQ(ierr);

    ierr = VecLoad(pGlobal, viewer); CHKERRQ(ierr);


    // destroy PetscViewer
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // write


PetscErrorCode Solutions::convert2Velocity(const Mat &RInv)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    Vec                 U;

    ierr = VecDuplicate(UGlobal, &U); CHKERRQ(ierr);

    ierr = MatMult(RInv, UGlobal, U); CHKERRQ(ierr);

    ierr = VecSwap(UGlobal, U); CHKERRQ(ierr);

    ierr = VecDestroy(&U); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


PetscErrorCode Solutions::convert2Flux(const Mat &R)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    Vec                 F;

    ierr = VecDuplicate(UGlobal, &F); CHKERRQ(ierr);

    ierr = MatMult(R, UGlobal, F); CHKERRQ(ierr);

    ierr = VecSwap(UGlobal, F); CHKERRQ(ierr);

    ierr = VecDestroy(&F); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


PetscErrorCode Solutions::applyIC(const FlowDescription &flow)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;
    std::vector<Vec>    UGlobalUnPacked(dim);

    ierr = DMCompositeGetAccessArray(mesh->UPack, UGlobal, dim,
            nullptr, UGlobalUnPacked.data()); CHKERRQ(ierr);

    for(PetscInt comp=0; comp<dim; ++comp)
    {
        ierr = VecSet(UGlobalUnPacked[comp], flow.IC[comp]); CHKERRQ(ierr);
    }

    ierr = DMCompositeRestoreAccessArray(mesh->UPack, UGlobal, dim,
            nullptr, UGlobalUnPacked.data()); CHKERRQ(ierr);

    ierr = VecSet(pGlobal, 0.0); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


std::ostream &operator<< (std::ostream &os, const Solutions &soln)
{
    os << soln.info;
    return os;
}
