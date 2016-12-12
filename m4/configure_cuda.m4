# CONFIGURE_CUDA65
# ---------------
# brief: Configures optional package CUDA required with AmgX.

AC_DEFUN([CONFIGURE_CUDA],[

echo
echo "================================="
echo "Configuring optional package CUDA"
echo "================================="

PACKAGE_SETUP_ENVIRONMENT

# check for presence of `--with-cuda=PATH`
AC_ARG_WITH([cuda],
            AS_HELP_STRING([--with-cuda=PATH],
                           [set CUDA directory]),
            [if test ! -d "$withval" ; then
               AC_MSG_ERROR([it is necessary to specify an existing directory 
                             when using --with-cuda=PATH])
             fi
             CUDA_DIR=$withval])

AC_MSG_NOTICE([using CUDA: $CUDA_DIR])

CUDA_INC_PATH=-I$CUDA_DIR/include
CUDA_LIB_PATH=-L$CUDA_DIR/lib64
CUDA_LIBRARY=-lcudart\ -lcublas\ -lcusparse

CPPFLAGS_PREPEND($CUDA_INC_PATH)
LDFLAGS_PREPEND($CUDA_LIB_PATH)
LIBS_PREPEND($CUDA_LIBRARY_PATH)

AC_CHECK_HEADER([cuda_runtime.h], ,
                AC_MSG_ERROR([could not find header file cuda_runtime.h]))

AC_MSG_CHECKING([for CUDA version])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[
#include <cuda_runtime.h>
#define CUDART_VERSION_EQ65 (CUDART_VERSION==6050)
]], [[
#if (CUDART_VERSION_EQ65)
#else
asdf
#endif
]]
                                  )],
                  [CUDART_VERSION_VALID=yes],
                  [CUDART_VERSION_VALID=no])
AC_MSG_RESULT([${CUDART_VERSION_VALID}])
if test "$CUDART_VERSION_VALID" = no; then
  AC_MSG_ERROR([invalid CUDA version detected; please use CUDA 6.5])
fi

AC_CHECK_LIB([cudart], [cudaRuntimeGetVersion], ,
             AC_MSG_ERROR([could not find library cudart]))
AC_CHECK_LIB([cublas], [cublasGetVersion], ,
             AC_MSG_ERROR([could not find library cublas]))
AC_CHECK_LIB([cusparse], [cusparseGetVersion], ,
             AC_MSG_ERROR([could not find library cusparse]))

PACKAGE_CPPFLAGS_PREPEND($CUDA_INC_PATH)
PACKAGE_LDFLAGS_PREPEND($CUDA_LIB_PATH)
PACKAGE_LIBS_PREPEND($CUDA_LIBRARY)

PACKAGE_RESTORE_ENVIRONMENT

AC_SUBST(CUDA_DIR, $CUDA_DIR)

echo

])
