/**
 * \file solutionsimple.cpp
 * \brief Implementation of the members of the class
 *        petibm::solution::SolutionSimple.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 */

#include <petibm/io.h>
#include <petibm/parser.h>
#include <petibm/solutionsimple.h>

namespace petibm
{
namespace solution
{
using namespace type;

// Default destructor.
SolutionSimple::~SolutionSimple() = default;

// Constructor; initialize the solution object.
SolutionSimple::SolutionSimple(const Mesh &inMesh)
{
    init(inMesh);
}  // SolutionSimple

// Create the flow field solutions.
PetscErrorCode SolutionSimple::init(const Mesh &inMesh)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    // make a shared pointer pointing to underlying mesh
    mesh = inMesh;

    // obtain MPI information from CartesianMesh object
    comm = mesh->comm;
    mpiSize = mesh->mpiSize;
    mpiRank = mesh->mpiRank;

    // get a copy of dim
    dim = mesh->dim;

    // create global composite vector for q
    ierr = DMCreateGlobalVector(mesh->UPack, &UGlobal); CHKERRQ(ierr);

    // create global composite vector for pressure
    ierr = DMCreateGlobalVector(mesh->da[3], &pGlobal); CHKERRQ(ierr);

    // create info string
    ierr = createInfoString(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // init

// Create a string with information about the solution.
PetscErrorCode SolutionSimple::createInfoString()
{
    PetscErrorCode ierr;
    std::stringstream ss;

    PetscFunctionBeginUser;

    if (mpiRank == 0)
    {
        PetscInt size;
        ss << std::string(80, '=') << std::endl;
        ss << "Solution Vectors:" << std::endl;
        ss << std::string(80, '=') << std::endl;
        ss << "\tDimension: " << dim << std::endl;
        ss << std::endl;
        ierr = VecGetSize(UGlobal, &size); CHKERRQ(ierr);
        ss << "\tLength of Global Packed Velocity Vector: " << size
           << std::endl;
        ss << std::endl;
        ierr = VecGetSize(pGlobal, &size); CHKERRQ(ierr);
        ss << "\tLength of Global Pressure Vector: " << size << std::endl;
        ss << std::endl;
    }

    info = ss.str();

    PetscFunctionReturn(0);
}  // createInfoString

// Convert velocity fluxes to velocity components.
PetscErrorCode SolutionSimple::convert2Velocity(const Mat &RInv)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    Vec U;  // temporary Vec for the velocity components
    ierr = VecDuplicate(UGlobal, &U); CHKERRQ(ierr);
    ierr = MatMult(RInv, UGlobal, U); CHKERRQ(ierr);
    ierr = VecSwap(UGlobal, U); CHKERRQ(ierr);
    ierr = VecDestroy(&U); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // convert2Velocity

// Convert velocity components to velocity fluxes.
PetscErrorCode SolutionSimple::convert2Flux(const Mat &R)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    Vec F;  // temporary Vec for the fluxes
    ierr = VecDuplicate(UGlobal, &F); CHKERRQ(ierr);
    ierr = MatMult(R, UGlobal, F); CHKERRQ(ierr);
    ierr = VecSwap(UGlobal, F); CHKERRQ(ierr);
    ierr = VecDestroy(&F); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // convert2Flux

// Set initial conditions of the flow fields.
PetscErrorCode SolutionSimple::setInitialConditions(const YAML::Node &node)
{
    PetscErrorCode ierr;
    type::RealVec1D ICs;
    std::vector<Vec> UGlobalUnpacked(dim);

    PetscFunctionBeginUser;

    // parse initial conditions for the velocity vector field
    ierr = parser::parseICs(node, ICs); CHKERRQ(ierr);

    // get individual velocity components from the packed Vec object
    ierr = DMCompositeGetAccessArray(mesh->UPack, UGlobal, dim, nullptr,
                                     UGlobalUnpacked.data()); CHKERRQ(ierr);

    for (PetscInt i = 0; i < dim; ++i)
    {
        ierr = VecSet(UGlobalUnpacked[i], ICs[i]); CHKERRQ(ierr);
    }

    ierr = DMCompositeRestoreAccessArray(mesh->UPack, UGlobal, dim, nullptr,
                                         UGlobalUnpacked.data());
    CHKERRQ(ierr);

    // pressure scalar field initialized with zeros
    ierr = VecSet(pGlobal, 0.0); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // setInitialConditions

// Write flow field solutions to a file.
PetscErrorCode SolutionSimple::write(const std::string &filePath) const
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    std::vector<Vec> vecs(dim + 1);
    std::vector<std::string> names(dim + 1);

    // set names for u, v, and/or w
    for (int f = 0; f < dim; ++f) names[f] = type::fd2str[type::Field(f)];
    // set name for pressure
    names.back() = "p";

    // get individual velocity components from the packed Vec object
    ierr = DMCompositeGetAccessArray(mesh->UPack, UGlobal, dim, nullptr,
                                     vecs.data()); CHKERRQ(ierr);
    // get a reference to pressure Vec
    vecs.back() = pGlobal;

    // write to a HDF5 file
    ierr = io::writeHDF5Vecs(comm, filePath, "/", names, vecs); CHKERRQ(ierr);

    // nullify the reference to pressure Vec
    vecs.back() = PETSC_NULL;
    // return individual Vec objects to the packed Vec object
    ierr = DMCompositeRestoreAccessArray(mesh->UPack, UGlobal, dim, nullptr,
                                         vecs.data()); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // write

// Read the flow field solutions from a file.
PetscErrorCode SolutionSimple::read(const std::string &filePath)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    std::vector<Vec> vecs(dim + 1);
    std::vector<std::string> names(dim + 1);

    // set names for u, v, and/or w
    for (int f = 0; f < dim; ++f) names[f] = type::fd2str[type::Field(f)];
    // set name for pressure
    names.back() = "p";

    // get individual velocity components from the packed Vec object
    ierr = DMCompositeGetAccessArray(mesh->UPack, UGlobal, dim, nullptr,
                                     vecs.data()); CHKERRQ(ierr);
    // get a reference to pressure Vec
    vecs.back() = pGlobal;

    // read Vec objects from a HDF5 file
    ierr = io::readHDF5Vecs(comm, filePath, "/", names, vecs); CHKERRQ(ierr);

    // nullify the reference to pressure Vec
    vecs.back() = PETSC_NULL;
    // return individual Vec objects to the packed Vec object
    ierr = DMCompositeRestoreAccessArray(mesh->UPack, UGlobal, dim, nullptr,
                                         vecs.data()); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}  // read

}  // end of namespace solution
}  // end of namespace petibm
