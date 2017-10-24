# CONFIGURE_CUDA
# ---------------
# brief: Configures optional package CUDA required with AmgX.

AC_DEFUN([CONFIGURE_CUDA],[

echo
echo "================================="
echo "Configuring optional package CUDA"
echo "================================="

PACKAGE_INITIALIZE_ENVIRONMENT
PACKAGE_SETUP_ENVIRONMENT

AC_ARG_VAR([CUDA_DIR], [CUDA directory])

AC_ARG_WITH([cuda-dir],
            AS_HELP_STRING([--with-cuda-dir=PATH],
                           [set CUDA directory]),
            [CUDA_DIR=$withval],
            [CUDA_DIR=/usr/local/cuda])

if test ! -d $CUDA_DIR ; then
  AC_MSG_ERROR([$CUDA_DIR is not a valid directory;
please use `--with-cuda-dir` to provide CUDA directory
or set the environment variable CUDA_DIR])
fi

AC_SUBST(CUDA_DIR, $CUDA_DIR)
AC_MSG_NOTICE([using CUDA: $CUDA_DIR])

CUDA_CPPFLAGS="-I$CUDA_DIR/include"
if test -d "$CUDA_DIR/lib"; then
  CUDA_LDFLAGS="-L$CUDA_DIR/lib -Wl,-rpath,$CUDA_DIR/lib"
elif test -d "$CUDA_DIR/lib64"; then
  CUDA_LDFLAGS="-L$CUDA_DIR/lib64 -Wl,-rpath,$CUDA_DIR/lib64"
fi
CUDA_LIBS="-lcudart -lcublas -lcusparse -lcusolver"
CPPFLAGS_PREPEND($CUDA_CPPFLAGS)
LDFLAGS_PREPEND($CUDA_LDFLAGS)

AC_CHECK_HEADER([cuda_runtime.h],
                [],
                [AC_MSG_ERROR([could not find header file cuda_runtime.h])])

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
  AC_MSG_ERROR([invalid CUDA version detected; please use CUDA 6.5 or above])
fi

AC_CHECK_LIB([cudart], [cudaRuntimeGetVersion], ,
             AC_MSG_ERROR([could not find library cudart; check config.log]))
AC_CHECK_LIB([cublas], [cublasGetVersion], ,
             AC_MSG_ERROR([could not find library cublas; check config.log]))
AC_CHECK_LIB([cusparse], [cusparseGetVersion], ,
             AC_MSG_ERROR([could not find library cusparse; check config.log]))

PACKAGE_RESTORE_ENVIRONMENT

echo

])
