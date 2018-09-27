/**
 * \file solutionsimple.h
 * \brief Definition of class petibm::solution::SolutionSimple.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 */

#pragma once

#include <petscsys.h>

#include <petibm/mesh.h>
#include <petibm/solution.h>

namespace petibm
{
namespace solution
{
/**
 * \brief Class to hold the velocity vector field and pressure scalar field.
 *
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

    // documentation: see petibm::solution::SolutionBase::setInitialConditions
    virtual PetscErrorCode setInitialConditions(const YAML::Node &node);

    // documentation: see petibm::solution::SolutionBase::convert2Velocity
    virtual PetscErrorCode convert2Velocity(const Mat &Rinv);

    // documentation: see petibm::solution::SolutionBase::convert2Flux
    virtual PetscErrorCode convert2Flux(const Mat &R);

    // documentation: see petibm::solution::SolutionBase::write
    virtual PetscErrorCode write(const std::string &filePath) const;

    // documentation: see petibm::solution::SolutionBase::read
    virtual PetscErrorCode read(const std::string &filePath);

protected:
    // documentation: see petibm::solution::SolutionBase::init
    virtual PetscErrorCode init(const type::Mesh &mesh);

    /**
     * \brief Create a string with information about the solution.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode createInfoString();

};  // SolutionSimple

}  // end of namespace solution
}  // end of namespace petibm
