/**
 * \file boundarysimple.h
 * \brief Definition of boundary::BoundarySimple
 * \author Anush Krishnan (anush@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \copyright MIT.
 */


# pragma once

// STL
# include <vector>

// PETSc
# include <petscvec.h>

// here goes headers from our PetIBM
# include <petibm/boundary.h>


namespace petibm
{
namespace boundary
{
/**
 * \brief An implementation of petibm::boundary::BoundaryBase.
 * \see boundaryModule, petibm::type::Boundary, petibm::boundary::BoundaryBase
 * \ingroup boundaryModule
 */
class BoundarySimple : public BoundaryBase
{
public:

    /** \copydoc BoundaryBase::BoundaryBase(const type::Mesh &, const YAML::Node &) */
    BoundarySimple(const type::Mesh &mesh, const YAML::Node &node);

    /** \copydoc BoundaryBase::~BoundaryBase */
    virtual ~BoundarySimple() = default;


    // implementation of BoundaryBase::setGhostICs
    virtual PetscErrorCode setGhostICs(const type::Solution &soln);

    // implementation of BoundaryBase::updateEqs
    virtual PetscErrorCode updateEqs(const type::Solution &soln, const PetscReal &dt);

    // implementation of BoundaryBase::updateGhostValues
    virtual PetscErrorCode updateGhostValues(const type::Solution &soln);

    // implementation of BoundaryBase::copyValues2LocalVecs
    virtual PetscErrorCode copyValues2LocalVecs(std::vector<Vec> &lclVecs) const;

protected:

    // implementation of BoundaryBase::init
    virtual PetscErrorCode init(const type::Mesh &mesh, const YAML::Node &node);

}; // BoundarySimple
} // end of namespace boundary
} // end of namespace petibm
