/***************************************************************************//**
 * \file solutions.h
 * \author Anush Krishnan (anush@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Definition of the class `SimulationParameters`.
 */


// PetIBM
# include <petibm/solutions.h>
# include <petibm/parser.h>


namespace petibm
{
namespace solution
{

using namespace type;

/** \copydoc Solutions::Solutions() */
Solutions::Solutions() = default;


/** \copydoc Solutions::~Solutions() */
Solutions::~Solutions() = default;


/** \copydoc Solutions::Solutions(const CartesianMesh &, const OutputType &) */
Solutions::Solutions(const Mesh &inMesh)
{
    init(inMesh);
}


PetscErrorCode Solutions::init(const Mesh &inMesh)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;

    // make a shared pointer pointing to associated mesh
    mesh = inMesh; 

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

    // create info string
    ierr = createInfoString(); CHKERRQ(ierr);

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


PetscErrorCode Solutions::applyIC(const YAML::Node &node)
{
    PetscFunctionBeginUser;

    PetscErrorCode      ierr;
    type::RealVec1D     values;
    std::vector<Vec>    UGlobalUnPacked(dim);
    
    ierr = parser::parseICs(node, values); CHKERRQ(ierr);

    ierr = DMCompositeGetAccessArray(mesh->UPack, UGlobal, dim,
            nullptr, UGlobalUnPacked.data()); CHKERRQ(ierr);

    for(PetscInt comp=0; comp<dim; ++comp)
    {
        ierr = VecSet(UGlobalUnPacked[comp], values[comp]); CHKERRQ(ierr);
    }

    ierr = DMCompositeRestoreAccessArray(mesh->UPack, UGlobal, dim,
            nullptr, UGlobalUnPacked.data()); CHKERRQ(ierr);

    ierr = VecSet(pGlobal, 0.0); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

} // end of namespace solution
} // end of namespace petibm
