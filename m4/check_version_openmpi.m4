# CHECKVERSION_OPENMPI
# ---------------
# brief: Check version of OpenMPI is compatible with AmgX.

AC_DEFUN([CHECK_VERSION_OPENMPI],[

echo
echo "========================"
echo "Check version of OpenMPI"
echo "========================"

PACKAGE_SETUP_ENVIRONMENT

AC_CHECK_HEADER([mpi.h], ,
                AC_MSG_ERROR([could not find header file mpi.h]))

AC_MSG_CHECKING([for OpenMPI])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[
#include <mpi.h>
]], [[
#ifdef OMPI_MPI_H
#else
asdf
#endif
]]
                                  )],
                  [IS_OMPI=yes],
                  [IS_OMPI=no])
AC_MSG_RESULT([${IS_OMPI}])
if test "$IS_OMPI" = no; then
  AC_MSG_ERROR([must use OpenMPI (AmgX works with OpenMPI)])
fi

AC_MSG_CHECKING([for version of OpenMPI])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[
#include <mpi.h>
# define CHECK (OMPI_MAJOR_VERSION >= 1 && OMPI_MINOR_VERSION >= 8)
]], [[
#if (CHECK)
#else
asdf
#endif
]]
                                  )],
                  [OMPI_VERSION_VALID=yes],
                  [OMPI_VERSION_VALID=no])
AC_MSG_RESULT([${OMPI_VERSION_VALID}])
if test "$OMPI_VERSION_VALID" = no; then
  AC_MSG_ERROR([OpenMPI should not be older than version 1.8 (AmgX currently only accept OpenMPI 1.8)])
fi

PACKAGE_RESTORE_ENVIRONMENT

echo

])
