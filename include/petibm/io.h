/**
 * \file io.h
 * \brief Prototypes of I/O functions.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
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
/** 
 * \brief A namespace of I/O-related functions.
 * \ingroup miscModule
 */
namespace io
{
    
/**
 * \brief Read the number and coordinates of Lagrangian points from a file.
 * \param file [in] the path to the file.
 * \param nPts [out] the total number of Lagrangian points.
 * \param coords [out] 2D STL vector of coordinates.
 * \return PetscErrorCode.
 * \ingroup miscModule
 */
PetscErrorCode readLagrangianPoints(
        const std::string &file, PetscInt &nPts, type::RealVec2D &coords);

/**
 * \brief Print information of "info" string of a parallel object to standard output.
 * \param info [in] a local string on this process.
 * \return PetscErrorCode.
 * \ingroup miscModule
 */
PetscErrorCode print(const std::string &info);


/**
 * \brief Write a vector of Vecs with names to a HDF5 file.
 * \param comm [in] MPI communicator (should be the same as the one in Vecs).
 * \param filePath [in] path of the file to write in.
 * \param loc [in] where in the HDF5 file that data will be written in.
 * \param names [in] a std::vector of std::string for the names of the Vecs.
 * \param vecs [in] a std::vector of Vecs.
 * \param mode [in] either FILE_MODE_WRITE (default) or FILE_MODE_APPEND.
 * \return PetscErrorCode.
 * \ingroup miscModule
 * 
 * Note: this function won't check the length of `vecs` and `names`.
 */
PetscErrorCode writeHDF5Vecs(const MPI_Comm comm, const std::string &filePath,
    const std::string &loc, const std::vector<std::string> &names,
    const std::vector<Vec> &vecs, const PetscFileMode mode=FILE_MODE_WRITE);


/**
 * \brief Write a vector of raw arrays with names to a HDF5 file.
 * \param comm [in] MPI communicator.
 * \param filePath [in] path of the file to write in.
 * \param loc [in] where in the HDF5 file that data will be written in.
 * \param names [in] a std::vector of std::string for the names of each array.
 * \param n [in] a std::vector of integer for the length of each array.
 * \param vecs [in] a std::vector of raw arrays (i.e., PetscReal*).
 * \param mode [in] either FILE_MODE_WRITE (default) or FILE_MODE_APPEND.
 * \return PetscErrorCode.
 * \ingroup miscModule
 * 
 * Note: this function won't check the length of `vecs` and `names`.
 */
PetscErrorCode writeHDF5Vecs(const MPI_Comm comm, const std::string &filePath,
    const std::string &loc, const std::vector<std::string> &names,
    const std::vector<PetscInt> &n, const std::vector<PetscReal*> &vecs,
    const PetscFileMode mode=FILE_MODE_WRITE);


/**
 * \brief Write type::RealVec2D with names to a HDF5 file.
 * \param comm [in] MPI communicator.
 * \param filePath [in] path of the file to write in.
 * \param loc [in] where in the HDF5 file that data will be written in.
 * \param names [in] a std::vector of std::string for the names of each RealVec1D.
 * \param vecs [in] a type::RealVec2D, i.e., std::vector<type::RealVec1D>.
 * \param mode [in] either FILE_MODE_WRITE (default) or FILE_MODE_APPEND.
 * \return PetscErrorCode.
 * \ingroup miscModule
 * 
 * Note: this function won't check the length of `vecs` and `names`.
 */
PetscErrorCode writeHDF5Vecs(const MPI_Comm comm, const std::string &filePath,
    const std::string &loc, const std::vector<std::string> &names,
    const type::RealVec2D &vecs, const PetscFileMode mode=FILE_MODE_WRITE);


/**
 * \brief read a vector of Vecs matching the provided names from a HDF5 file.
 * \param comm [in] MPI communicator (should be the same as the one in Vecs).
 * \param filePath [in] path of the file to read from.
 * \param loc [in] where in the HDF5 file that data will be read from.
 * \param names [in] a std::vector of std::string for the names of the Vecs.
 * \param vecs [out] a std::vector of Vecs.
 * \return PetscErrorCode.
 * \ingroup miscModule
 * 
 * Note: this function won't check the length of `vecs` and `names`.
 */
PetscErrorCode readHDF5Vecs(const MPI_Comm comm, const std::string &filePath,
    const std::string &loc, const std::vector<std::string> &names,
    std::vector<Vec> &vecs);


/**
 * \brief Write a summary of the PETSc logging into a ASCII file.
 * \param comm [in] MPI communicator.
 * \param filePath [in] Path of the file to write in.
 * \return PetscErrorCode
 * \ingroup miscModule
 */
PetscErrorCode writePetscLog(const MPI_Comm comm, const std::string &filePath);

} // end of io
} // end of petibm
