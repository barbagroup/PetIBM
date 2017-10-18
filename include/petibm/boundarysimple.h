/***************************************************************************//**
 * \file boundarysimple.h
 * \author Anush Krishnan (anus@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Definition of the class `BoundarySimple`.
 */


# pragma once

// STL
# include <vector>

// PETSc
# include <petscvec.h>

// here goes headers from our PetIBM
# include <petibm/singleboundary.h>
# include <petibm/boundary.h>


namespace petibm
{
namespace boundary
{
/** \brief a class holding data of ghosted boundaries. */
class BoundarySimple : public BoundaryBase
{
public:

    /** \copydoc petibm::boundary::BoundaryBase::BoundaryBase */
    BoundarySimple(const type::Mesh &mesh, const YAML::Node &node);

    /** \copydoc petibm::boundary::BoundaryBase::~BoundaryBase */
    virtual ~BoundarySimple();


    /** \copydoc petibm::boundary::BoundaryBase::setGhostICs */
    virtual PetscErrorCode setGhostICs(const type::Solution &soln);

    /** \copydoc petibm::boundary::BoundaryBase::updateEqs */
    virtual PetscErrorCode updateEqs(const type::Solution &soln, const PetscReal &dt);

    /** \copydoc petibm::boundary::BoundaryBase::updateGhostValues */
    virtual PetscErrorCode updateGhostValues(const type::Solution &soln);

    /** \copydoc petibm::boundary::BoundaryBase::copyValues2LocalVecs */
    virtual PetscErrorCode copyValues2LocalVecs(std::vector<Vec> &lclVecs) const;

protected:

    /** \copydoc petibm::boundary::BoundaryBase::init */
    virtual PetscErrorCode init(const type::Mesh &mesh, const YAML::Node &node);

    /** \brief a 2D vector holding all single boundaries. */
    std::vector<std::vector<type::SingleBoundary>>  bds;

};
} // end of namespace boundary
} // end of namespace petibm
