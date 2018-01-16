# CONFIGURE_AMGX
# ---------------
# brief: Configures required package AmgX.

AC_DEFUN([CONFIGURE_AMGX],[

CONFIGURE_CUDA

echo
echo "================================="
echo "Configuring optional package AmgX"
echo "================================="

PACKAGE_INITIALIZE_ENVIRONMENT
PACKAGE_SETUP_ENVIRONMENT

AC_ARG_VAR([AMGX_DIR], [AmgX directory])

AC_ARG_WITH([amgx-dir],
            AS_HELP_STRING([--with-amgx-dir=PATH],
                           [set AmgX directory]),
            [AMGX_DIR=$withval],
            [])

if test ! -d $AMGX_DIR || test "x$AMGX_DIR" = "x" ; then
  AC_MSG_ERROR([$AMGX_DIR is not a valid AmgX directory;
please use `--with-amgx-dir` to provide AmgX directory
or set the environment variable AMGX_DIR;
otherwise, disable amgx])
fi

AC_SUBST(AMGX_DIR, $AMGX_DIR)
AC_MSG_NOTICE([using AmgX: $AMGX_DIR])

if test -f "$AMGX_DIR/base/include/amgx_c.h"; then
  AMGX_CPPFLAGS="-DHAVE_AMGX -I$AMGX_DIR/base/include $CUDA_CPPFLAGS"
elif test -f "$AMGX_DIR/include/amgx_c.h"; then
  AMGX_CPPFLAGS="-DHAVE_AMGX -I$AMGX_DIR/include $CUDA_CPPFLAGS"
elif test -f "$AMGX_DIR/examples/amgx_c.h"; then
  AMGX_CPPFLAGS="-DHAVE_AMGX -I$AMGX_DIR/examples $CUDA_CPPFLAGS"
fi
if test -d "$AMGX_DIR/lib"; then
  AMGX_LDFLAGS="-L$AMGX_DIR/lib -Wl,-rpath,$AMGX_DIR/lib $CUDA_LDFLAGS"
elif test -d "$AMGX_DIR/lib64"; then
  AMGX_LDFLAGS="-L$AMGX_DIR/lib64 -Wl,-rpath,$AMGX_DIR/lib64 $CUDA_LDFLAGS"
fi
CPPFLAGS_PREPEND($AMGX_CPPFLAGS)
LDFLAGS_PREPEND($AMGX_LDFLAGS)

# check for presence of header file amgx_c.h
AC_CHECK_HEADER([amgx_c.h], ,
                AC_MSG_ERROR([could not find header file amgx_c.h]))

if test "x$enable_shared" = "xyes"; then
  AMGX_LIBS="-lamgxsh $CUDA_LIBS"
  # check for presence of library
  AC_CHECK_LIB([amgxsh],
               [AMGX_initialize], ,
               AC_MSG_ERROR([could not use library amgxsh; check config.log]))
else
  AMGX_LIBS="-lamgx $CUDA_LIBS"
  # check for presence of library
  AC_CHECK_LIB([amgx],
               [AMGX_initialize], ,
               AC_MSG_ERROR([could not use library amgx; check config.log]),
               [$CUDA_LIBS -lcusolver -fopenmp])
fi

PACKAGE_RESTORE_ENVIRONMENT

echo

CONFIGURE_AMGXWRAPPER

]) # CONFIGURE_AMGX
