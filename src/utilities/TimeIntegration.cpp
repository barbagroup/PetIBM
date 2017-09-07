/***************************************************************************//**
 * \file SimulationParameters.h
 * \author Anush Krishnan (anush@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Definition of members of the class `TimeIntegration`.
 */


// here goes headers from our PetIBM
#include "TimeIntegration.h"


namespace petibm
{
namespace utilities
{

using namespace types;

/** \copydoc TimeIntegration::TimeIntegration(). */
TimeIntegration::TimeIntegration() = default;


/** \copydoc TimeIntegration::~TimeIntegration(). */
TimeIntegration::~TimeIntegration() = default;


/** \copydoc TimeIntegration::TimeIntegration(const types::TimeScheme &). */
TimeIntegration::TimeIntegration(const TimeScheme &method) { init(method); }


/** \copydoc TimeIntegration::init(const types::TimeScheme &). */
PetscErrorCode TimeIntegration::init(const TimeScheme &method)
{
    PetscFunctionBeginUser;

    // hard copy the scheme
    scheme = method;


    // set up coefficients.
    switch (scheme)
    {
        case NONE:
            implicitCoeff = 0.0;
            nExplicit = 0;
            break;
        case EULER_EXPLICIT:
            implicitCoeff = 0.0; // n+1 coefficient
            nExplicit = 1;
            explicitCoeffs.push_back(1.0); // n coefficient
            break;
        case EULER_IMPLICIT:
            implicitCoeff = 1.0; // n+1 coefficient
            nExplicit = 0;
            break;
        case ADAMS_BASHFORTH_2:
            implicitCoeff = 0.0;  // n+1 coefficient
            nExplicit = 2;
            explicitCoeffs.push_back(1.5);  // n coefficient
            explicitCoeffs.push_back(-0.5); // n-1 coefficient
            break;
        case CRANK_NICOLSON:
            implicitCoeff = 0.5; // n+1 coefficient
            nExplicit = 1;
            explicitCoeffs.push_back(0.5); // n coefficient
            break;
    }

    PetscFunctionReturn(0);
}

} // end of namespace utilities
} // end of namespace petibm
