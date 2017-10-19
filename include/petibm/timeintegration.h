/***************************************************************************//**
 * \file SimulationParameters.h
 * \author Anush Krishnan (anush@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Definition of TimeIntegration related code.
 */


# pragma once

// STL
# include <string>
# include <memory>

// PETSc headers
# include <petscsys.h>

// YAML CPP
# include <yaml-cpp/yaml.h>

// headers from PetIBM
# include <petibm/type.h>


namespace petibm
{
namespace timeintegration
{

/**
* \class TimeIntegrationBase
* \brief Stores information about temporal integration schemes.
*/
class TimeIntegrationBase
{

public:
    
    /** \brief name of current instance. */
    const std::string                 name;

    /** \brief name of the scheme. */
    const std::string                 scheme;

    /** \brief coefficient of inplicit term. */
    const PetscReal                   implicitCoeff;

    /** \brief number of explicit terms. */
    const PetscInt                    nExplicit;

    /** \brief coefficients of explicit terms. */
    const type::RealVec1D             explicitCoeffs;


    /** \brief default constructor. */
    TimeIntegrationBase(): TimeIntegrationBase("none", "none", 0.0, 0, {}) {};
    
    /**
     * \brief constructor (normally not being used publicly).
     *
     * \param inName [in] the name of the instance.
     * \param inScheme [in] the name of the scheme.
     * \param inImplicitCoeff [in] implicit coefficient.
     * \param inNEcplicit [in] number of explicit coefficients.
     * \param inExplicitCoeffs [in] a std::vector holding all explit coefficients.
     */
    TimeIntegrationBase(
            const std::string &inName,
            const std::string &inScheme,
            const PetscReal &inImplicitCoeff,
            const PetscInt  &inNEcplicit,
            const type::RealVec1D &inExplicitCoeffs):
        name(inName), scheme(inScheme), implicitCoeff(inImplicitCoeff),
        nExplicit(inNEcplicit), explicitCoeffs(inExplicitCoeffs) {};

    /** \brief destructor. */
    virtual ~TimeIntegrationBase() = default;

}; // TimeIntegrationBase


/** \brief 1st order explicit Euler. */
class Euler_Explicit : public TimeIntegrationBase
{
    Euler_Explicit(const std::string &name):
        TimeIntegrationBase(name, "1st order explicit Euler", 0.0, 1, {1.0}) {};
    
    ~Euler_Explicit() = default;
};


/** \brief 1st order implicit Euler. */
class Euler_Implicit : public TimeIntegrationBase
{
    Euler_Implicit(const std::string &name):
        TimeIntegrationBase(name, "1st order implicit Euler", 1.0, 0, {}) {};
    
    ~Euler_Implicit() = default;
};


/** \brief 2nd order Adams-Bashforth. */
class Adams_Bashforth_2 : public TimeIntegrationBase
{
    Adams_Bashforth_2(const std::string &name):
        TimeIntegrationBase(name, "2nd order Adams-Bashforth", 0.0, 2, {1.5, -0.5}) {};
    
    ~Adams_Bashforth_2() = default;
};


/** \brief 2nd order Crank-Nicolson. */
class Crank_Nicolson : public TimeIntegrationBase
{
    Crank_Nicolson(const std::string &name):
        TimeIntegrationBase(name, "2nd order Crank-Nicolson", 0.5, 1, {0.5}) {};
    
    ~Crank_Nicolson() = default;
};

} // end of namespace timeintegration


namespace type
{
    /** \brief definition of type::TimeIntegration. */
    typedef std::shared_ptr<timeintegration::TimeIntegrationBase> TimeIntegration;
}


namespace timeintegration
{
    /**
     * \brief factory function for type::TimeIntegration.
     *
     * \param name [in] name of the instance.
     * \param node [in] YAML::Node of all configuration.
     * \param integration [out] resulting TimeIntegration object.
     *
     * \return PetscErrorCode.
     */
    PetscErrorCode createTimeIntegration(
            const std::string &name, const YAML::Node &node,
            type::TimeIntegration &integration);
}

} // end of namespace petibm
