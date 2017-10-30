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


/**
 * \brief wirte a vector of Vecs with names to a HDF5 file.
 *
 * \param comm [in] MPI communicator (should be the same as the one in Vecs).
 * \param file [in] file path + file name (without extension).
 * \param names [in] a std::vector of std::string for the names of the Vecs.
 * \param vecs [in] a std::vector of Vecs.
 * 
 * Note: this function won't check the length of `vecs` and `names`.
 * 
 * \param mode [in] either FILE_MODE_WRITE (default) or FILE_MODE_APPEND.
 *
 * \return PetscErrorCode.
 */
PetscErrorCode writeHDF5Vecs(const MPI_Comm comm, const std::string &file, 
        const std::vector<std::string> &names, const std::vector<Vec> &vecs, 
        const PetscFileMode mode=FILE_MODE_WRITE);


/**
 * \brief read a vector of Vecs matching the provided names from a HDF5 file.
 *
 * \param comm [in] MPI communicator (should be the same as the one in Vecs).
 * \param file [in] file path + file name (without extension).
 * \param names [in] a std::vector of std::string for the names of the Vecs.
 * \param vecs [out] a std::vector of Vecs.
 * 
 * Note: this function won't check the length of `vecs` and `names`.
 * 
 * \return PetscErrorCode.
 */
PetscErrorCode readHDF5Vecs(const MPI_Comm comm, const std::string &file,
        const std::vector<std::string> &names, std::vector<Vec> &vecs);

} // end of io
} // end of petibm
