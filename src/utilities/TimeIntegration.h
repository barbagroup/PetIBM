/***************************************************************************//**
 * \file SimulationParameters.h
 * \author Anush Krishnan (anush@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Definition of the class `TimeIntegration`.
 */


# pragma once

// here goes C++ STL
# include <vector>

// here goes PETSc headers
# include <petscsys.h>

// here goes headers from our PetIBM
# include "types.h"


namespace petibm
{
namespace utilities
{

/**
* \class TimeIntegration
* \brief Stores information about temporal integration schemes.
*/
class TimeIntegration
{

public:

    /** \brief type of time-stepping scheme. */
    types::TimeScheme           scheme;

    /** \brief coefficient of inplicit term. */
    PetscReal                   implicitCoeff;

    /** \brief number of explicit terms. */
    PetscInt                    nExplicit;

    /** \brief coefficients of explicit terms. */
    types::RealVec1D            explicitCoeffs;


    /** \brief default constructor. */
    TimeIntegration();

    /** \brief constructor. */
    TimeIntegration(const types::TimeScheme &method);

    /** \brief destructor. */
    ~TimeIntegration();

    /** \brief initialization. */
    PetscErrorCode init(const types::TimeScheme &method);

}; // TimeIntegration

} // end of namespace utilities
} // end of namespace petibm
