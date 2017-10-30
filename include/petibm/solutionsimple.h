/***************************************************************************//**
 * \file solutions.h
 * \author Anush Krishnan (anush@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Definition of the class `Solutions`.
 */


# pragma once

// PETSc
# include <petscsys.h>

// PetIBM
# include <petibm/mesh.h>
# include <petibm/solution.h>


namespace petibm
{
namespace solution
{
/** \brief a class holding only velocity and pressure solution. */
class SolutionSimple : public SolutionBase
{
public:

    /** \copydoc petibm::solution::SolutionBase::SolutionBase */
    SolutionSimple(const type::Mesh &mesh);

    /** \copydoc petibm::solution::SolutionBase::~SolutionBase */
    virtual ~SolutionSimple();

    /** \copydoc petibm::solution::SolutionBase::applyIC */
    virtual PetscErrorCode applyIC(const YAML::Node &node);

    /** \copydoc petibm::solution::SolutionBase::convert2Velocity */
    virtual PetscErrorCode convert2Velocity(const Mat &Rinv);

    /** \copydoc petibm::solution::SolutionBase::convert2Flux */
    virtual PetscErrorCode convert2Flux(const Mat &R);
    
    /** \copydoc petibm::solution::SolutionBase::write */
    virtual PetscErrorCode write(const std::string &file) const;
    
    /** \copydoc petibm::solution::SolutionBase::read */
    virtual PetscErrorCode read(const std::string &file);


protected:

    /** \copydoc petibm::solution::SolutionBase::init */
    virtual PetscErrorCode init(const type::Mesh &mesh);

    /**
     * \brief create a string for information.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode createInfoString();
};
} // end of namespace solution
} // end of namespace petibm
