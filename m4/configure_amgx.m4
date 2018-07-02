# CONFIGURE_AMGX
# --------------
# brief: Configures required package AmgX (required with AmgXWrapper).

AC_DEFUN([CONFIGURE_AMGX], [

CONFIGURE_CUDA

echo
echo "================================="
echo "Configuring optional package AmgX"
echo "================================="

PACKAGE_INITIALIZE_ENVIRONMENT
PACKAGE_SETUP_ENVIRONMENT

AC_ARG_WITH([amgx-dir],
            AS_HELP_STRING([--with-amgx-dir=PATH],
                           [set AmgX directory]),
            [AMGX_DIR=$withval],
            [])

AC_ARG_WITH([amgx-include],
            AS_HELP_STRING([--with-amgx-include=PATH],
                           [set AmgX include directory]),
            [AMGX_INCLUDE_DIR=$withval],
            [])

AC_ARG_WITH([amgx-lib],
            AS_HELP_STRING([--with-amgx-lib=PATH],
                           [set AmgX lib directory]),
            [AMGX_LIB_DIR=$withval],
            [])


if test "x$AMGX_DIR" != "x"; then
  CHECK_AMGX_DIR($AMGX_DIR)
else
  if test "x$AMGX_INCLUDE_DIR" != "x"; then
    CHECK_AMGX_INCLUDE_DIR($AMGX_INCLUDE_DIR)
  fi
  if test "x$AMGX_LIB_DIR" != "x"; then
    CHECK_AMGX_LIB_DIR($AMGX_LIB_DIR)
  fi
fi

AC_CHECK_HEADER([amgx_c.h], [],
                AC_MSG_ERROR([could not find amgx_c.h; check config.log]))

if test "x$enable_shared" = "xyes"; then
  AMGX_LIBS="-lamgxsh $CUDA_LIBS"
  AC_CHECK_LIB([amgxsh], [AMGX_initialize], [],
               AC_MSG_ERROR([could not use library amgxsh; check config.log]))
elif test "x$enable_static" = "xyes"; then
  AMGX_LIBS="-lamgx $CUDA_LIBS"
  AC_CHECK_LIB([amgx], [AMGX_initialize], [],
               AC_MSG_ERROR([could not use library amgx; check config.log]),
               [$CUDA_LIBS])
fi

PACKAGE_RESTORE_ENVIRONMENT

echo

]) # CONFIGURE_AMGX


# CHECK_AMGX_DIR
# --------------
# brief: Checks the existence of the AmgX directory.
AC_DEFUN([CHECK_AMGX_DIR], [

AMGX_DIR=$1
AC_MSG_NOTICE([using AmgX: $AMGX_DIR])

if test ! -d "$AMGX_DIR"; then
  AC_MSG_ERROR([$AMGX_DIR is not a valid directory...
Use '--with-amgx-dir=PATH' to provide the directory of the package.])
fi
CHECK_AMGX_INCLUDE_DIR($AMGX_DIR/include)
if test -d "$AMGX_DIR/lib64"; then
  CHECK_AMGX_LIB_DIR($AMGX_DIR/lib64)
else
  CHECK_AMGX_LIB_DIR($AMGX_DIR/lib)
fi

]) # CHECK_AMGX_DIR


# CHECK_AMGX_INCLUDE_DIR
# ----------------------
# brief: Checks the existence of the AmgX include directory.
AC_DEFUN([CHECK_AMGX_INCLUDE_DIR], [

AMGX_INCLUDE_DIR=$1
AC_MSG_NOTICE([using AmgX include directory: $AMGX_INCLUDE_DIR])

if test -d "$AMGX_INCLUDE_DIR"; then
  AMGX_CPPFLAGS="-DHAVE_AMGX -I$AMGX_INCLUDE_DIR $CUDA_CPPFLAGS"
  AC_SUBST(AMGX_CPPFLAGS, $AMGX_CPPFLAGS)
  CPPFLAGS_APPEND($AMGX_CPPFLAGS)
else
  AC_MSG_ERROR([$AMGX_INCLUDE_DIR is not a valid directory...
Use '--with-amgx-include=PATH' to provide the include directory of the package.])
fi

]) # CHECK_AMGX_INCLUDE_DIR


# CHECK_AMGX_LIB_DIR
# ------------------
# brief: Checks the existence of the AmgX lib directory.
AC_DEFUN([CHECK_AMGX_LIB_DIR], [

AMGX_LIB_DIR=$1
AC_MSG_NOTICE([using AmgX lib directory: $AMGX_LIB_DIR])

if test -d "$AMGX_LIB_DIR"; then
  AMGX_LDFLAGS="-L$AMGX_LIB_DIR -Wl,-rpath,$AMGX_LIB_DIR $CUDA_LDFLAGS"
  AC_SUBST(AMGX_LDFLAGS, $AMGX_LDFLAGS)
  LDFLAGS_APPEND($AMGX_LDFLAGS)
else
  AC_MSG_ERROR([$AMGX_LIB_DIR is not a valid directory...
Use '--with-amgx-lib=PATH' to provide the lib directory of the package.])
fi

]) # CHECK_AMGX_LIB_DIR
