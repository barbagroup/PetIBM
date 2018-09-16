/**
 * \file timeintegration.h
 * \brief Definition of TimeIntegration related classes.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 */

#pragma once

// STL
#include <memory>
#include <string>

// PETSc headers
#include <petscsys.h>

// YAML CPP
#include <yaml-cpp/yaml.h>

// headers from PetIBM
#include <petibm/type.h>

/**
 * \defgroup timeModule Time-integration schemes
 * \brief Objects holding informations of time-integration schemes.
 *
 * API users should use petibm::timeintegration::createTimeIntegration to
 * create desired instances.
 *
 * \see petibm::timeintegration::createTimeIntegration,
 * petibm::type::TimeIntegration \ingroup petibm
 */

namespace petibm
{
namespace timeintegration
{
/**
 * \class TimeIntegrationBase
 * \brief Base (abstract) class that stores information of temporal integration.
 * \see timeModule
 * \ingroup timeModule
 */
class TimeIntegrationBase
{
public:
    /** \brief Name of current instance. */
    const std::string name;

    /** \brief Name of the scheme. */
    const std::string scheme;

    /** \brief Coefficient of implicit term. */
    const PetscReal implicitCoeff;

    /** \brief Number of explicit terms. */
    const PetscInt nExplicit;

    /** \brief Coefficients of explicit terms. */
    const type::RealVec1D explicitCoeffs;

    /** \brief Default constructor. */
    TimeIntegrationBase() : TimeIntegrationBase("none", "none", 0.0, 0, {}){};

    /**
     * \brief Constructor (normally not being used publicly).
     * \param inName [in] the name of the instance.
     * \param inScheme [in] the name of the scheme.
     * \param inImplicitCoeff [in] implicit coefficient.
     * \param inNEcplicit [in] number of explicit coefficients.
     * \param inExplicitCoeffs [in] a std::vector holding all explicit
     * coefficients.
     */
    TimeIntegrationBase(const std::string &inName, const std::string &inScheme,
                        const PetscReal &inImplicitCoeff,
                        const PetscInt &inNEcplicit,
                        const type::RealVec1D &inExplicitCoeffs)
        : name(inName),
          scheme(inScheme),
          implicitCoeff(inImplicitCoeff),
          nExplicit(inNEcplicit),
          explicitCoeffs(inExplicitCoeffs){};

    /** \brief Destructor. */
    virtual ~TimeIntegrationBase() = default;

    /**
     * \brief Print information to standard output.
     * \return PetscErrorCode.
     */
    PetscErrorCode printInfo() const;

};  // TimeIntegrationBase

/**
 * \brief An implementation of TimeIntegrationBase for 1st order explicit Euler.
 * \see timeModule
 * \ingroup timeModule
 */
class Euler_Explicit : public TimeIntegrationBase
{
public:
    /**
     * \brief Constructor.
     * \param name [in] the name of the instance.
     */
    Euler_Explicit(const std::string &name)
        : TimeIntegrationBase(name, "1st order explicit Euler", 0.0, 1,
                              {1.0}){};

    /** \copydoc TimeIntegrationBase::~TimeIntegrationBase */
    virtual ~Euler_Explicit() = default;
};  // Euler_Explicit

/**
 * \brief An implementation of TimeIntegrationBase for 1st order implicit Euler.
 * \see timeModule
 * \ingroup timeModule
 */
class Euler_Implicit : public TimeIntegrationBase
{
public:
    /**
     * \brief Constructor.
     * \param name [in] the name of the instance.
     */
    Euler_Implicit(const std::string &name)
        : TimeIntegrationBase(name, "1st order implicit Euler", 1.0, 0, {}){};

    /** \copydoc TimeIntegrationBase::~TimeIntegrationBase */
    virtual ~Euler_Implicit() = default;
};  // Euler_Implicit

/**
 * \brief An implementation of TimeIntegrationBase for 2nd order
 * Adams-Bashforth. \see timeModule \ingroup timeModule
 */
class Adams_Bashforth_2 : public TimeIntegrationBase
{
public:
    /**
     * \brief Constructor.
     * \param name [in] the name of the instance.
     */
    Adams_Bashforth_2(const std::string &name)
        : TimeIntegrationBase(name, "2nd order Adams-Bashforth", 0.0, 2,
                              {1.5, -0.5}){};

    /** \copydoc TimeIntegrationBase::~TimeIntegrationBase */
    virtual ~Adams_Bashforth_2() = default;
};  // Adams_Bashforth_2

/**
 * \brief An implementation of TimeIntegrationBase for 2nd order Crank-Nicolson.
 * \see timeModule
 * \ingroup timeModule
 */
class Crank_Nicolson : public TimeIntegrationBase
{
public:
    /**
     * \brief Constructor.
     * \param name [in] the name of the instance.
     */
    Crank_Nicolson(const std::string &name)
        : TimeIntegrationBase(name, "2nd order Crank-Nicolson", 0.5, 1,
                              {0.5}){};

    /** \copydoc TimeIntegrationBase::~TimeIntegrationBase */
    virtual ~Crank_Nicolson() = default;
};  // Crank_Nicolson

}  // end of namespace timeintegration

namespace type
{
/**
 * \brief Definition of type::TimeIntegration.
 * \see timeModule, petibm::timeintegration::createTimeIntegration
 * \ingroup timeModule
 */
typedef std::shared_ptr<timeintegration::TimeIntegrationBase> TimeIntegration;
}  // end of namespace type

namespace timeintegration
{
/**
 * \brief factory function for type::TimeIntegration.
 * \param name [in] name of the instance.
 * \param node [in] YAML::Node of all configuration.
 * \param integration [out] resulting TimeIntegration object.
 * \return PetscErrorCode.
 * \see timeModule, petibm::type::TimeIntegration
 * \ingroup timeModule
 */
PetscErrorCode createTimeIntegration(const std::string &name,
                                     const YAML::Node &node,
                                     type::TimeIntegration &integration);
}  // end of namespace timeintegration

}  // end of namespace petibm
