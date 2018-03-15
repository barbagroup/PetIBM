# CONFIGURE_AMGXWRAPPER
# ---------------
# brief: Configures optional package AmgXWrapper.

AC_DEFUN([CONFIGURE_AMGXWRAPPER],[

echo
echo "========================================"
echo "Configuring required package AmgXWrapper"
echo "========================================"

PACKAGE_INITIALIZE_ENVIRONMENT
PACKAGE_SETUP_ENVIRONMENT

AC_ARG_VAR([AMGXWRAPPER_DIR], [AmgXWrapper directory])

AC_ARG_WITH([amgxwrapper-dir],
            AS_HELP_STRING([--with-amgxwrapper-dir=PATH],
                           [set AmgXWrapper directory]),
            [AMGXWRAPPER_DIR=$withval],
            [])

AC_ARG_ENABLE([amgxwrapper],
              AS_HELP_STRING([--enable-amgxwrapper],
                             [download and install AmgXWrapper-1.3]),
              [IS_INSTALL_AMGXWRAPPER=yes],
              [IS_INSTALL_AMGXWRAPPER=no])

if test "x$AMGXWRAPPER_DIR" != "x"; then
  if test -d "$AMGXWRAPPER_DIR"; then
    SETUP_AMGXWRAPPER
  else
    AC_MSG_ERROR([
$AMGXWRAPPER_DIR is not a valid directory;
please use '--with-amgxwrapper-dir=PATH' to provide the directory of the package])
  fi
else
  if test "x$IS_INSTALL_AMGXWRAPPER" = "xyes"; then
    INSTALL_AMGXWRAPPER
  fi
fi

CPPFLAGS_APPEND($AMGX_CPPFLAGS)
CPPFLAGS_APPEND($PETSC_CPPFLAGS)
AC_CHECK_HEADER([AmgXSolver.hpp],
                [],
                [AC_MSG_ERROR([Couldn't find AmgXSolver.hpp...
Please use '--with-amgxwrapper-dir=PATH' to provide the directory of the package
or '--enable-amgxwrapper' to download and install AmgXWrapper-1.3 (AmgX-2.0).
])])

AMGXWRAPPER_LIBS="-lAmgXWrapper $AMGX_LIBS"
AC_SUBST(AMGXWRAPPER_LIBS, $AMGXWRAPPER_LIBS)

PACKAGE_RESTORE_ENVIRONMENT

echo

])


AC_DEFUN([SETUP_AMGXWRAPPER], [

AC_SUBST(AMGXWRAPPER_DIR, $AMGXWRAPPER_DIR)
AC_MSG_NOTICE([using AmgXWrapper: $AMGXWRAPPER_DIR])
AMGXWRAPPER_CPPFLAGS="-I$AMGXWRAPPER_DIR/include $AMGX_CPPFLAGS"
AC_SUBST(AMGXWRAPPER_CPPFLAGS, $AMGXWRAPPER_CPPFLAGS)
if test -d "$AMGXWRAPPER_DIR/lib"; then
  AMGXWRAPPER_LDFLAGS="-L$AMGXWRAPPER_DIR/lib -Wl,-rpath,$AMGXWRAPPER_DIR/lib $AMGX_LDFLAGS"
elif test -d "$AMGXWRAPPER_DIR/lib64"; then
  AMGXWRAPPER_LDFLAGS="-L$AMGXWRAPPER_DIR/lib64 -Wl,-rpath,$AMGXWRAPPER_DIR/lib64 $AMGX_LDFLAGS"
fi
AC_SUBST(AMGXWRAPPER_LDFLAGS, $AMGXWRAPPER_LDFLAGS)
CPPFLAGS_PREPEND($AMGXWRAPPER_CPPFLAGS)
LDFLAGS_PREPEND($AMGXWRAPPER_LDFLAGS)

]) # SETUP_AMGXWRAPPER


AC_DEFUN([INSTALL_AMGXWRAPPER], [

PACKAGE_INITIALIZE_ENVIRONMENT
PACKAGE_SETUP_ENVIRONMENT
AC_MSG_NOTICE([

 ****************************************
 *      Downloading and installing      *
 *           AmgXWrapper-1.4            *
 ****************************************
])
echo "*** INFO *** Downloading AmgXWrapper-1.4... "
VERSION=1.4
TARBALL=v$VERSION.tar.gz
URL=https://github.com/barbagroup/AmgXWrapper/archive/$TARBALL
wget $URL -P /tmp
AMGXWRAPPER_DIR=$BUILDDIR/externalpackages/AmgXWrapper-$VERSION
AC_SUBST(AMGXWRAPPER_DIR, $AMGXWRAPPER_DIR)
mkdir -p $AMGXWRAPPER_DIR/build
tar -xzf /tmp/$TARBALL -C $AMGXWRAPPER_DIR --strip-components=1
rm -f /tmp/$TARBALL
echo "*** INFO *** Building AmgXWrapper-1.4... "
cd $AMGXWRAPPER_DIR/build
if test "x$enable_shared" = "xyes"; then
  cmake $AMGXWRAPPER_DIR \
    -DAMGX_DIR=$AMGX_DIR \
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
    -DAMGX_DIR=$AMGX_DIR \
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
AC_SUBST(AMGXWRAPPER_CPPFLAGS, $AMGX_CPPFLAGS)
AC_SUBST(AMGXWRAPPER_LDFLAGS, $AMGX_LDFLAGS)
PACKAGE_RESTORE_ENVIRONMENT

echo

]) # INSTALL_AMGXWRAPPER
