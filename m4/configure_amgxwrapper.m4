# CONFIGURE_AMGXWRAPPER
# ---------------------
# brief: Configures optional package AmgXWrapper.
AC_DEFUN([CONFIGURE_AMGXWRAPPER], [

CONFIGURE_AMGX

echo
echo "========================================"
echo "Configuring required package AmgXWrapper"
echo "========================================"

PACKAGE_INITIALIZE_ENVIRONMENT
PACKAGE_SETUP_ENVIRONMENT

AC_ARG_WITH([amgxwrapper-dir],
            AS_HELP_STRING([--with-amgxwrapper-dir=PATH],
                           [set AmgXWrapper directory]),
            [AMGXWRAPPER_DIR=$withval],
            [])

AC_ARG_WITH([amgxwrapper-include],
            AS_HELP_STRING([--with-amgxwrapper-include=PATH],
                           [set AmgXWrapper include directory]),
            [AMGXWRAPPER_INCLUDE_DIR=$withval],
            [])

AC_ARG_WITH([amgxwrapper-lib],
            AS_HELP_STRING([--with-amgxwrapper-lib=PATH],
                           [set AmgXWrapper lib directory]),
            [AMGXWRAPPER_LIB_DIR=$withval],
            [])

AC_ARG_ENABLE([amgxwrapper],
              AS_HELP_STRING([--enable-amgxwrapper],
                             [download and install AmgXWrapper-1.4]),
              [IS_INSTALL_AMGXWRAPPER=yes],
              [IS_INSTALL_AMGXWRAPPER=no])

if test "x$AMGXWRAPPER_DIR" != "x"; then
  CHECK_AMGXWRAPPER_DIR($AMGXWRAPPER_DIR)
else
  if test "x$AMGXWRAPPER_INCLUDE_DIR" != "x"; then
    CHECK_AMGXWRAPPER_INCLUDE_DIR($AMGXWRAPPER_INCLUDE_DIR)
  fi
  if test "x$AMGXWRAPPER_LIB_DIR" != "x"; then
    CHECK_AMGXWRAPPER_LIB_DIR($AMGXWRAPPER_LIB_DIR)
  fi
  if test "x$IS_INSTALL_AMGXWRAPPER" = "xyes"; then
    INSTALL_AMGXWRAPPER
  fi
fi

AC_CHECK_HEADER([AmgXSolver.hpp], [],
                [AC_MSG_ERROR([could not find AmgXSolver.hpp...
Use '--with-amgxwrapper-dir=PATH' to provide the directory of the package or
Use '--with-amgxwrapper-include=PATH' and '--with-amgxwrapper-lib=PATH' or
Use '--enable-amgxwrapper' to download, build, and install AmgXWrapper-1.4.])])

AMGXWRAPPER_LIBS="-lAmgXWrapper $AMGX_LIBS"
AC_SUBST(AMGXWRAPPER_LIBS, $AMGXWRAPPER_LIBS)

PACKAGE_RESTORE_ENVIRONMENT

echo

])


# CHECK_AMGXWRAPPER_DIR
# ---------------------
# brief: Checks the existence of the AmgXWrapper directory.
AC_DEFUN([CHECK_AMGXWRAPPER_DIR], [

AMGXWRAPPER_DIR=$1
AC_MSG_NOTICE([using AmgXWrapper: $AMGXWRAPPER_DIR])

if test ! -d "$AMGXWRAPPER_DIR"; then
  AC_MSG_ERROR([$AMGXWRAPPER_DIR is not a valid directory...
Use '--with-amgxwrapper-dir=PATH' to provide the directory of the package.])
fi
CHECK_AMGXWRAPPER_INCLUDE_DIR($AMGXWRAPPER_DIR/include)
if test -d "$AMGXWRAPPER_DIR/lib64"; then
  CHECK_AMGXWRAPPER_LIB_DIR($AMGXWRAPPER_DIR/lib64)
else
  CHECK_AMGXWRAPPER_LIB_DIR($AMGXWRAPPER_DIR/lib)
fi

]) # CHECK_AMGXWRAPPER_DIR


# CHECK_AMGX_INCLUDE_DIR
# ----------------------
# brief: Checks the existence of the AmgXWrapper include directory.
AC_DEFUN([CHECK_AMGXWRAPPER_INCLUDE_DIR], [

AMGXWRAPPER_INCLUDE_DIR=$1
AC_MSG_NOTICE([using AmgXWrapper include directory: $AMGXWRAPPER_INCLUDE_DIR])

if test -d "$AMGXWRAPPER_INCLUDE_DIR"; then
  AMGXWRAPPER_CPPFLAGS="-I$AMGXWRAPPER_INCLUDE_DIR $AMGX_CPPFLAGS $PETSC_CPPFLAGS"
  AC_SUBST(AMGXWRAPPER_CPPFLAGS, $AMGXWRAPPER_CPPFLAGS)
  CPPFLAGS_APPEND($AMGXWRAPPER_CPPFLAGS)
else
  AC_MSG_ERROR([$AMGXWRAPPER_INCLUDE_DIR is not a valid directory...
Use '--with-amgxwrapper-include=PATH' to provide the include directory of the package.])
fi

]) # CHECK_AMGXWRAPPER_INCLUDE_DIR


# CHECK_AMGXWRAPPER_LIB_DIR
# -------------------------
# brief: Checks the existence of the AmgXWrapper lib directory.
AC_DEFUN([CHECK_AMGXWRAPPER_LIB_DIR], [

AMGXWRAPPER_LIB_DIR=$1
AC_MSG_NOTICE([using AmgXWrapper lib directory: $AMGXWRAPPER_LIB_DIR])

if test -d "$AMGXWRAPPER_LIB_DIR"; then
  AMGXWRAPPER_LDFLAGS="-L$AMGXWRAPPER_LIB_DIR -Wl,-rpath,$AMGXWRAPPER_LIB_DIR $AMGX_LDFLAGS $PETSC_LDFLAGS"
  AC_SUBST(AMGXWRAPPER_LDFLAGS, $AMGXWRAPPER_LDFLAGS)
  LDFLAGS_APPEND($AMGXWRAPPER_LDFLAGS)
else
  AC_MSG_ERROR([$AMGXWRAPPER_LIB_DIR is not a valid directory...
Use '--with-amgxwrapper-lib=PATH' to provide the lib directory of the package.])
fi

]) # CHECK_AMGXWRAPPER_LIB_DIR


# INSTALL_AMGXWRAPPER
# -------------------------
# brief: Downloads, builds, and installs AmgXWrapper.
AC_DEFUN([INSTALL_AMGXWRAPPER], [

version=1.4
tarball=v$version.tar.gz
url=https://github.com/barbagroup/AmgXWrapper/archive/$tarball
AMGXWRAPPER_DIR=$BUILDDIR/externalpackages/AmgXWrapper-$version
if test ! -d "$AMGXWRAPPER_DIR"; then
  echo "downloading AmgXWrapper-$version... "
  wget -q $url -P /tmp
  mkdir -p $AMGXWRAPPER_DIR/build
  tar -xzf /tmp/$tarball -C $AMGXWRAPPER_DIR --strip-components=1
  rm -f /tmp/$tarball
fi
echo "building and installing AmgXWrapper-$version... "
cd $AMGXWRAPPER_DIR/build
if test "x$enable_shared" = "xyes"; then
  cmake $AMGXWRAPPER_DIR \
    -DAMGX_INCLUDE_DIRS=$AMGX_INCLUDE_DIR \
    -DAMGX_LIBRARIES=$AMGX_LIB_DIR/libamgxsh.so \
    -DCUDA_DIR=$CUDA_DIR \
    -DPETSC_DIR=$PETSC_DIR \
    -DPETSC_ARCH=$PETSC_ARCH \
    -DBUILD_SHARED_LIBS=ON \
    -DCMAKE_INSTALL_PREFIX=$prefix \
    -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON
  make all
  make install
  rm -f $AMGXWRAPPER_DIR/build/CMakeCache.txt
fi
if test "x$enable_static" = "xyes"; then
  cmake $AMGXWRAPPER_DIR \
    -DAMGX_INCLUDE_DIRS=$AMGX_INCLUDE_DIR \
    -DAMGX_LIBRARIES=$AMGX_LIB_DIR/libamgx.a \
    -DCUDA_DIR=$CUDA_DIR \
    -DPETSC_DIR=$PETSC_DIR \
    -DPETSC_ARCH=$PETSC_ARCH \
    -DBUILD_SHARED_LIBS=OFF \
    -DCMAKE_INSTALL_PREFIX=$prefix \
    -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON
  make all
  make install
  rm -f $AMGXWRAPPER_DIR/build/CMakeCache.txt
fi
cd $BUILDDIR
AMGXWRAPPER_CPPFLAGS="$AMGX_CPPFLAGS $PETSC_CPPFLAGS"
AC_SUBST(AMGXWRAPPER_CPPFLAGS, $AMGXWRAPPER_CPPFLAGS)
CPPFLAGS_APPEND($AMGXWRAPPER_CPPFLAGS)
AMGXWRAPPER_LDFLAGS="$AMGX_LDFLAGS $PETSC_LDFLAGS"
AC_SUBST(AMGXWRAPPER_LDFLAGS, $AMGXWRAPPER_LDFLAGS)
LDFLAGS_APPEND($AMGXWRAPPER_LDFLAGS)

echo

]) # INSTALL_AMGXWRAPPER
