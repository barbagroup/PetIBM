/*
 * io.h
 * Copyright (C) 2017 Pi-Yueh Chuang <pychuang@gwu.edu>
 *
 * Distributed under terms of the MIT license.
 */

# pragma once

// STL
# include <string>

// PETSc
# include <petscsys.h>

// PetIBM
# include <petibm/type.h>
# include <petibm/singlebody.h>


namespace petibm
{
namespace io
{
    
/**
 * \brief read number and coordinates of Lagrangian points from a file.
 *
 * \param file [in] the path to the file.
 * \param nPts [out] the total number of Lagragian points.
 * \param coords [out] 2D STL vector of coordinates.
 *
 * \return PetscErrorCode.
 */
PetscErrorCode readLagrangianPoints(
        const std::string &file, PetscInt &nPts, type::RealVec2D &coords);

/**
 * \brief print information of "info" string of a parallel object to standard output.
 *
 * \param info [in] a local string on this process.
 *
 * \return PetscErrorCode.
 */
PetscErrorCode print(const std::string &info);

} // end of io
} // end of petibm
