# file: configure.ac
# author: Olivier Mesnard (mesnardo@gwu.edu)
# brief: Input file for Autoconf.

# set-up
AC_PREREQ([2.69])
AC_INIT([PetIBM], [0.2.0dev], [mesnardo@gwu.edu, anush@bu.edu])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])
AC_CONFIG_SRCDIR([src/PetIBM.cpp])
AC_CONFIG_HEADERS([config/config.h])
AM_SILENT_RULES([yes])
AM_MAINTAINER_MODE([disable])
AC_LANG([C++])

# check for programs
AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXX
AC_PROG_CXXCPP
AC_PROG_AWK
AC_PROG_INSTALL
AC_PROG_RANLIB
AM_PROG_AR
AC_PROG_MKDIR_P
AC_PROG_CC_C_O
CONFIGURE_DOXYGEN

# check for typedefs, structures, and compiler characteristics
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_INT64_T
AC_CHECK_TYPES([ptrdiff_t])

# check for header files
AC_CHECK_HEADERS([arpa/inet.h fcntl.h float.h limits.h netdb.h stddef.h stdlib.h string.h strings.h sys/param.h sys/time.h sys/timeb.h unistd.h wchar.h wctype.h])

# check for library functions
AC_FUNC_FORK
AC_FUNC_MALLOC
AC_FUNC_MMAP
AC_FUNC_STRERROR_R
AC_CHECK_FUNCS([dup2 fchdir getcwd getpagesize gettimeofday memset mkdir munmap pow regcomp rmdir socket sqrt strcasecmp strchr strdup strerror strrchr strstr strtol strtoull])

# check for optional and required third-party libraries
PACKAGE_INITIALIZE_ENVIRONMENT
CONFIGURE_PETSC
AC_ARG_WITH([amgx],, [CHECK_VERSION_OPENMPI
                      CONFIGURE_CUDA
                      CONFIGURE_AMGX])
AM_CONDITIONAL([WITH_AMGXWRAPPER], [test "x$WITH_AMGX" = xyes])
PACKAGE_SETUP_ENVIRONMENT

# list of Makefiles to generate
AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/utilities/Makefile
                 src/solvers/Makefile
                 external/Makefile
                 external/yaml-cpp-0.5.1/Makefile
                 external/gtest-1.7.0/Makefile
                 external/AmgXWrapper-1.0-beta/Makefile
                 tests/Makefile
                 tests/CartesianMesh/Makefile
                 tests/NavierStokes/Makefile
                 tests/TairaColonius/Makefile
                 tests/convectiveTerm/Makefile
                 tests/diffusiveTerm/Makefile
                 examples/Makefile])

# output message
AC_OUTPUT
echo
echo "====================="
echo "Configuration summary"
echo "====================="
AC_MSG_NOTICE([PetIBM appears to have been configured successfully

Prefix: '${prefix}'

CXX: ${CXX}
CXXFLAGS: ${CXXFLAGS}
LDFLAGS: ${LDFLAGS}

LIBS: ${LIBS}

PETSC_DIR: ${PETSC_DIR}
PETSC_ARCH: ${PETSC_ARCH}

Now type 'make @<:@<target>@:>@'
  where the optional @<:@<target>@:>@ is:
    all            - build all binaries
    install        - install everything
])
