/**
 * \file solutionsimple.h
 * \brief Definition of class solution::SolutionSimple.
 * \author Anush Krishnan (anus@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \copyright MIT.
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
/**
 * \brief A class holding only velocity and pressure solutions.
 * \see solutionModule, petibm::type::Solution, petibm::solution::createSolution
 * \ingroup solutionModule
 */
class SolutionSimple : public SolutionBase
{
public:

    /** \copydoc SolutionBase(const type::Mesh &) */
    SolutionSimple(const type::Mesh &mesh);

    /** \copydoc ~SolutionBase */
    virtual ~SolutionSimple();

    // doc is the same as solution::SolutionBase::applyIC
    virtual PetscErrorCode applyIC(const YAML::Node &node);

    // doc is the same as solution::SolutionBase::convert2Velocity
    virtual PetscErrorCode convert2Velocity(const Mat &Rinv);

    // doc is the same as solution::SolutionBase::convert2Flux
    virtual PetscErrorCode convert2Flux(const Mat &R);
    
    // doc is the same as solution::SolutionBase::write
    virtual PetscErrorCode write(const std::string &file) const;
    
    // doc is the same as solution::SolutionBase::read
    virtual PetscErrorCode read(const std::string &file);


protected:

    // doc is the same as SolutionBase::init */
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
