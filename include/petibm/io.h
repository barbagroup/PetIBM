/**
 * \file io.h
 * \brief Prototypes of I/O functions.
 * \copyright Copyright (c) 2016-2018, Barba group. All rights reserved.
 * \license BSD 3-Clause License.
 */

#pragma once

#include <string>

#include <petscsys.h>

#include <petibm/singlebody.h>
#include <petibm/type.h>

namespace petibm
{
/**
 * \brief A namespace of I/O-related functions.
 * \ingroup miscModule
 */
namespace io
{
/**
 * \brief Read the number Lagrangian points and their coordinates from a file.
 *
 * \param file [in] Path of the file to read from.
 * \param nPts [out] Number of Lagrangian points.
 * \param coords [out] Coordinates of the Lagrangian points as a 2D STL vector.
 *
 * \ingroup miscModule
 */
PetscErrorCode readLagrangianPoints(const std::string &file, PetscInt &nPts,
                                    type::RealVec2D &coords);

/**
 * \brief Print information of a parallel object to standard output.
 *
 * \param info [in] a local string on this process.
 *
 * \ingroup miscModule
 */
PetscErrorCode print(const std::string &info);

/**
 * \brief Write a vector of Vec objects to a HDF5 file.
 *
 * \param comm [in] MPI communicator (should be the same as the one in Vecs).
 * \param filePath [in] Path of the file to write in.
 * \param loc [in] Location in the HDF5 file for data to write.
 * \param names [in] Vector with the name of each Vec object to write.
 * \param vecs [in] Vector of Vec objects to write.
 * \param mode [in] Either FILE_MODE_WRITE (default) or FILE_MODE_APPEND.
 *
 * \ingroup miscModule
 *
 * Note: this function doesn't check the length of the vector of Vec objects
 * and of the vector of names.
 */
PetscErrorCode writeHDF5Vecs(const MPI_Comm comm, const std::string &filePath,
                             const std::string &loc,
                             const std::vector<std::string> &names,
                             const std::vector<Vec> &vecs,
                             const PetscFileMode mode = FILE_MODE_WRITE);

/**
 * \brief Write a vector of raw arrays to a HDF5 file.
 *
 * \param comm [in] MPI communicator.
 * \param filePath [in] Path of the file to write in.
 * \param loc [in] Location in the HDF5 file of data to write.
 * \param names [in] Vector with the name of each raw array to write.
 * \param n [in] Vector with the length of each raw array to write.
 * \param vecs [in] Vector of raw arrays to write.
 * \param mode [in] Either FILE_MODE_WRITE (default) or FILE_MODE_APPEND.
 *
 * \ingroup miscModule
 *
 * Note: this function doesn't check the length of the vector of raw arrays
 * and of the vector of names.
 */
PetscErrorCode writeHDF5Vecs(const MPI_Comm comm, const std::string &filePath,
                             const std::string &loc,
                             const std::vector<std::string> &names,
                             const std::vector<PetscInt> &n,
                             const std::vector<PetscReal *> &vecs,
                             const PetscFileMode mode = FILE_MODE_WRITE);

/**
 * \brief Write petibm::type::RealVec2D to a HDF5 file.
 *
 * \param comm [in] MPI communicator.
 * \param filePath [in] Path of the file to write in.
 * \param loc [in] Location in the HDF5 file of the data to write.
 * \param names [in] Vector with the name of each RealVec1D to write.
 * \param vecs [in] RealVec2D (i.e., std::vector<type::RealVec1D>) to write.
 * \param mode [in] Either FILE_MODE_WRITE (default) or FILE_MODE_APPEND.
 *
 * \ingroup miscModule
 *
 * Note: this function doesn't check the length of the RealVec2D object
 * and of the vector of names.
 */
PetscErrorCode writeHDF5Vecs(const MPI_Comm comm, const std::string &filePath,
                             const std::string &loc,
                             const std::vector<std::string> &names,
                             const type::RealVec2D &vecs,
                             const PetscFileMode mode = FILE_MODE_WRITE);

/**
 * \brief Read a vector of Vec objects from a HDF5 file.
 *
 * \param comm [in] MPI communicator (should be the same as the one in Vecs).
 * \param filePath [in] Path of the file to read from.
 * \param loc [in] Location in the HDF5 file of the data to read.
 * \param names [in] Vector with the name of each Vec object to read.
 * \param vecs [out] Vector of Vec objects to fill.
 *
 * \ingroup miscModule
 *
 * Note: this function doesn't check the length of the vector of Vec objects
 * and of the vector of names.
 */
PetscErrorCode readHDF5Vecs(const MPI_Comm comm, const std::string &filePath,
                            const std::string &loc,
                            const std::vector<std::string> &names,
                            std::vector<Vec> &vecs);

/**
 * \brief Write a summary of the PETSc logging into a ASCII file.
 *
 * \param comm [in] MPI communicator.
 * \param filePath [in] Path of the file to write in.
 *
 * \ingroup miscModule
 */
PetscErrorCode writePetscLog(const MPI_Comm comm, const std::string &filePath);

}  // namespace io

}  // namespace petibm
