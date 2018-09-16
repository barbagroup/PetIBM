/**
 * \file solutionsimple.h
 * \brief Definition of class solution::SolutionSimple.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 */

#pragma once

// PETSc
#include <petscsys.h>

// PetIBM
#include <petibm/mesh.h>
#include <petibm/solution.h>

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
    virtual PetscErrorCode write(const std::string &filePath) const;

    // doc is the same as solution::SolutionBase::read
    virtual PetscErrorCode read(const std::string &filePath);

protected:
    // doc is the same as SolutionBase::init */
    virtual PetscErrorCode init(const type::Mesh &mesh);

    /**
     * \brief create a string for information.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode createInfoString();
};  // SolutionSimple
}  // end of namespace solution
}  // end of namespace petibm
