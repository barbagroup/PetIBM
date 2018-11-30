# CONFIGURE_CUDA
# --------------
# brief: Configures optional package CUDA (required with AmgX).
AC_DEFUN([CONFIGURE_CUDA], [

echo
echo "================================="
echo "Configuring optional package CUDA"
echo "================================="

PACKAGE_INITIALIZE_ENVIRONMENT
PACKAGE_SETUP_ENVIRONMENT

AC_ARG_WITH([cuda-dir],
            AS_HELP_STRING([--with-cuda-dir=PATH],
                           [set CUDA directory]),
            [CUDA_DIR=$withval],
            [])

if test "x$CUDA_DIR" != "x"; then
  CHECK_CUDA_DIR($CUDA_DIR)
fi

AC_CHECK_HEADER([cuda_runtime.h], [],
                [AC_MSG_ERROR([could not find cuda_runtime.h; check config.log])])

AC_MSG_CHECKING([for CUDA version])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[
#include <cuda_runtime.h>
#define CUDART_VERSION_GE65 (CUDART_VERSION>=6050)
]], [[
#if (CUDART_VERSION_GE65)
#else
asdf
#endif
]]
                                  )],
                  [CUDART_VERSION_VALID=yes],
                  [CUDART_VERSION_VALID=no])
AC_MSG_RESULT([${CUDART_VERSION_VALID}])
if test "$CUDART_VERSION_VALID" = no; then
  AC_MSG_ERROR([invalid CUDA version detected; use CUDA-6.5 or above])
fi

CUDA_LIBS="-lcudart -lcublas -lcusparse -lcusolver"
AC_CHECK_LIB([cudart], [cudaRuntimeGetVersion], [],
             [AC_MSG_ERROR([could not find library cudart; check config.log])])
AC_CHECK_LIB([cublas], [cublasGetVersion], [],
             [AC_MSG_ERROR([could not find library cublas; check config.log])])
AC_CHECK_LIB([cusparse], [cusparseGetVersion], [],
             [AC_MSG_ERROR([could not find library cusparse; check config.log])])

PACKAGE_RESTORE_ENVIRONMENT

echo

]) # CONFIGURE_CUDA


# CHECK_CUDA_DIR
# --------------
# brief: Checks the existence of the CUDA directory.
AC_DEFUN([CHECK_CUDA_DIR], [

CUDA_DIR=$1
AC_MSG_NOTICE([using CUDA: $CUDA_DIR])

if test ! -d "$CUDA_DIR"; then
  AC_MSG_ERROR([$CUDA_DIR is not a valid directory...
Use '--with-cuda-dir=PATH' to provide the directory of the package.])
fi
CHECK_CUDA_INCLUDE_DIR($CUDA_DIR/include)
if test -d "$CUDA_DIR/lib64"; then
  CHECK_CUDA_LIB_DIR($CUDA_DIR/lib64)
else
  CHECK_CUDA_LIB_DIR($CUDA_DIR/lib)
fi

]) # CHECK_CUDA_DIR


# CHECK_CUDA_INCLUDE_DIR
# ----------------------
# brief: Checks the existence of the CUDA include directory.
AC_DEFUN([CHECK_CUDA_INCLUDE_DIR], [

CUDA_INCLUDE_DIR=$1
AC_MSG_NOTICE([using CUDA include directory: $CUDA_INCLUDE_DIR])

if test -d "$CUDA_INCLUDE_DIR"; then
  CUDA_CPPFLAGS="-I$CUDA_INCLUDE_DIR"
  AC_SUBST(CUDA_CPPFLAGS, $CUDA_CPPFLAGS)
  CPPFLAGS_APPEND($CUDA_CPPFLAGS)
else
  AC_MSG_ERROR([$CUDA_INCLUDE_DIR is not a valid directory...
Use '--with-cuda-include=PATH' to provide the include directory of the package.])
fi

]) # CHECK_CUDA_INCLUDE_DIR


# CHECK_CUDA_LIB_DIR
# ------------------
# brief: Checks the existence of the CUDA lib directory.
AC_DEFUN([CHECK_CUDA_LIB_DIR], [

CUDA_LIB_DIR=$1
AC_MSG_NOTICE([using CUDA lib directory: $CUDA_LIB_DIR])

if test -d "$CUDA_LIB_DIR"; then
  CUDA_LDFLAGS="-L$CUDA_LIB_DIR -Wl,-rpath,$CUDA_LIB_DIR"
  AC_SUBST(CUDA_LDFLAGS, $CUDA_LDFLAGS)
  LDFLAGS_APPEND($CUDA_LDFLAGS)
else
  AC_MSG_ERROR([$CUDA_LIB_DIR is not a valid directory...
Use '--with-cuda-lib=PATH' to provide the lib directory of the package.])
fi

]) # CHECK_CUDA_LIB_DIR
