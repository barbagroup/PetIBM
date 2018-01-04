# CONFIGURE_PETSC
# ---------------
# brief: Configures required package PETSc.

AC_DEFUN([CONFIGURE_PETSC], [

echo
echo "=================================="
echo "Configuring required package PETSc"
echo "=================================="

PACKAGE_INITIALIZE_ENVIRONMENT
PACKAGE_SETUP_ENVIRONMENT

# check for presence of `--with-petsc-dir=PATH` and `--with-petsc-arch=PATH`
AC_ARG_WITH([petsc-dir],
            AS_HELP_STRING([--with-petsc-dir=PATH],
                           [set PETSc directory]),
            [PETSC_DIR=$withval],
            [])
AC_ARG_WITH([petsc-arch],
            AS_HELP_STRING([--with-petsc-arch=PATH],
                           [set PETSc arch folder]),
            [PETSC_ARCH=$withval],
            [])

if test ! -d "$PETSC_DIR" ; then
  AC_MSG_ERROR([$PETSC_DIR is not a valid directory;
please use `--with-petsc-dir` to provide a PETSc directory])
else
  if test ! -d "$PETSC_DIR/$PETSC_ARCH" ; then
    AC_MSG_ERROR([$PETSC_DIR/$PETSC_ARCH is not a valid directory;
please use `--with-petsc-arch` to provide a PETSc arch])
  fi
fi

AC_MSG_NOTICE([using PETSc: $PETSC_DIR])
AC_MSG_NOTICE([with arch: $PETSC_ARCH])

# check for presence of file petscvariables
PETSCVARIABLES="$PETSC_DIR/$PETSC_ARCH/lib/petsc/conf/petscvariables"
AC_CHECK_FILE([$PETSCVARIABLES],
              [HAVE_PETSCVARIABLES=yes],
              [HAVE_PETSCVARIABLES=no])
if test "$HAVE_PETSCVARIABLES" = no; then
  AC_MSG_ERROR([could not find file petscvariables;
please use PETSc-3.7])
fi

PETSC_CC_INCLUDES=`grep "PETSC_CC_INCLUDES =" $PETSCVARIABLES | sed -e 's/.*=//' -e 's/^[ \t]*//'`
PETSC_EXTERNAL_LIB_BASIC=`grep "PETSC_EXTERNAL_LIB_BASIC =" $PETSCVARIABLES | sed -e 's/.*=//' -e 's/^[ \t]*//'`
PETSC_WITH_EXTERNAL_LIB=`grep "PETSC_WITH_EXTERNAL_LIB =" $PETSCVARIABLES | sed -e 's/.*=//' -e 's/^[ \t]*//'`

CPPFLAGS_PREPEND($PETSC_CC_INCLUDES)
AC_SUBST(PETSC_CPPFLAGS, $PETSC_CC_INCLUDES)
LIBS_PREPEND($PETSC_WITH_EXTERNAL_LIB)
AC_SUBST(PETSC_LDFLAGS, $PETSC_WITH_EXTERNAL_LIB)

# check for presence of header file petsc.h
AC_CHECK_HEADER([petsc.h], , 
                AC_MSG_ERROR([could not find header file petsc.h]))

# check for PETSc version
AC_MSG_CHECKING([for PETSc version])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[
#include <petscversion.h>

#ifndef PETSC_VERSION_GE
#define PETSC_VERSION_GE(MAJOR, MINOR, SUBMINOR) (!PETSC_VERSION_LT(MAJOR, MINOR, SUBMINOR))
#endif
]], [[
#if ((PETSC_VERSION_GE(3, 8, 0) && PETSC_VERSION_LT(3, 8, 4)) || !PETSC_VERSION_RELEASE)
#else
asdf
#endif
]]
                                  )],
                  [PETSC_VERSION_VALID=yes],
                  [PETSC_VERSION_VALID=no])
AC_MSG_RESULT([${PETSC_VERSION_VALID}])
if test "$PETSC_VERSION_VALID" = no; then
  AC_MSG_ERROR([invalid PETSc version detected; please use PETSc 3.8])
fi

PACKAGE_RESTORE_ENVIRONMENT

echo

]) # CONFIGURE_PETSC
