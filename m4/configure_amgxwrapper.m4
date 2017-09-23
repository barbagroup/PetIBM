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
                             [download and install AmgXWrapper-1.2]),
              [DOWNLOAD_AMGXWRAPPER=yes],
              [DOWNLOAD_AMGXWRAPPER=no])

if test "x$AMGXWRAPPER_DIR" != "x"; then
  if test -d "$AMGXWRAPPER_DIR"; then
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
  else
    AC_MSG_ERROR([
$AMGXWRAPPER_DIR is not a valid directory;
please use '--with-amgxwrapper-dir=PATH' to provide the directory of the package])
  fi
else
  if test "x$DOWNLOAD_AMGXWRAPPER" = "xyes"; then
    AC_MSG_NOTICE([

 ****************************************
 *      Downloading and installing      *
 *           AmgXWrapper-1.2            *
 ****************************************
])
    PACKAGE_INITIALIZE_ENVIRONMENT
    PACKAGE_SETUP_ENVIRONMENT
    echo "downloading AmgXWrapper-1.2... "
    VERSION=1.2
    TARBALL=v$VERSION.tar.gz
    URL=https://github.com/barbagroup/AmgXWrapper/archive/$TARBALL
    wget $URL -P /tmp
    AMGXWRAPPER_DIR=$BUILDDIR/externalpackages/AmgXWrapper-$VERSION
    AC_SUBST(AMGXWRAPPER_DIR, $AMGXWRAPPER_DIR)
    mkdir -p $AMGXWRAPPER_DIR/build
    tar -xzf /tmp/$TARBALL -C $AMGXWRAPPER_DIR --strip-components=1
    rm -f /tmp/$TARBALL
    echo "building AmgXWrapper-1.2... "
    cd $AMGXWRAPPER_DIR/build
    for VAL in "ON" "OFF"
    do
      cmake $AMGXWRAPPER_DIR \
        -DAMGX_DIR=$AMGX_DIR \
        -DCUDA_DIR=$CUDA_DIR \
        -DPETSC_DIR=$PETSC_DIR \
        -DPETSC_ARCH=$PETSC_ARCH \
        -DBUILD_SHARED_LIBS=$VAL \
        -DCMAKE_INSTALL_PREFIX=$prefix \
        -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON
      make all
      make install
    done
    echo "done! "
    cd $BUILDDIR
    AC_SUBST(AMGXWRAPPER_CPPFLAGS, $AMGX_CPPFLAGS)
    AC_SUBST(AMGXWRAPPER_LDFLAGS, $AMGX_LDFLAGS)
    PACKAGE_RESTORE_ENVIRONMENT
  fi
fi

CPPFLAGS_APPEND($AMGX_CPPFLAGS)
CPPFLAGS_APPEND($PETSC_CPPFLAGS)

AC_CHECK_HEADER([AmgXSolver.hpp],
                [],
                [AC_MSG_ERROR([Couldn't find AmgXSolver.hpp...
Please use '--with-amgxwrapper-dir=PATH' to provide the directory of the package
or '--enable-amgxwrapper' to download and install AmgXWrapper-1.1 (AmgX-2.0).
])])

AMGXWRAPPER_LIBS="-lAmgXWrapper"
AC_SUBST(AMGXWRAPPER_LIBS, $AMGXWRAPPER_LIBS)

PACKAGE_RESTORE_ENVIRONMENT

echo

])
