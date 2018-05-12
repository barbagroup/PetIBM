/**
 * \file io.hpp
 * \brief prototypes of functions in io.cpp
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \date 2017-05-31
 */


# pragma once


// PETSc
# include <petscmat.h>
# include <petscvec.h>


// StructArgs
# include "StructArgs.hpp"


/**
 * \brief print a header containing case info.
 *
 * \param args [in] a StructArgs instance.
 *
 * \return PetscErrorCode.
 */
PetscErrorCode printHeader(StructArgs &args);


/**
 * \brief print a footer.
 *
 * \param args [in] a StructArgs instance.
 *
 * \return PetscErrorCode.
 */
PetscErrorCode printFooter(StructArgs &args);
