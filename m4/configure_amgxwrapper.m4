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
                             [download and install AmgXWrapper-1.1 (AmgX-2.0)]),
              [DOWNLOAD_AMGXWRAPPER=yes],
              [DOWNLOAD_AMGXWRAPPER=no])

if test "x$AMGXWRAPPER_DIR" != "x"; then
  if test -d "$AMGXWRAPPER_DIR"; then
    AC_SUBST(AMGXWRAPPER_DIR, $AMGXWRAPPER_DIR)
    AC_MSG_NOTICE([using AmgXWrapper: $AMGXWRAPPER_DIR])
    AMGXWRAPPER_CPPFLAGS="-I$AMGXWRAPPER_DIR/include $AMGX_CPPFLAGS"
    AC_SUBST(AMGXWRAPPER_CPPFLAGS, $AMGXWRAPPER_CPPFLAGS)
    AMGXWRAPPER_LDFLAGS="-L$AMGXWRAPPER_DIR/lib -Wl,-rpath,$AMGXWRAPPER_DIR/lib $AMGX_LDFLAGS"
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
 *      AmgXWrapper-1.1 (AmgX-2.0)      *
 ****************************************
])
    PACKAGE_INITIALIZE_ENVIRONMENT
    PACKAGE_SETUP_ENVIRONMENT
    echo "downloading AmgXWrapper-1.1 (AmgX-2.0)... "
    ZIP_FILE=amgx-2.0.zip
    URL=https://github.com/mesnardo/AmgXWrapper/archive/$ZIP_FILE
    wget $URL -P /tmp
    unzip /tmp/$ZIP_FILE -d $BUILDDIR/externalpackages
    rm -f /tmp/$ZIP_FILE
    echo "building AmgXWrapper-1.1 (AmgX-2.0)... "
    AMGXWRAPPER_DIR=$BUILDDIR/externalpackages/AmgXWrapper-amgx-2.0
    AC_SUBST(AMGXWRAPPER_DIR, $AMGXWRAPPER_DIR)
    mkdir -p $AMGXWRAPPER_DIR/build
    cd $AMGXWRAPPER_DIR/build
    for VAL in "ON" "OFF"
    do
      cmake $AMGXWRAPPER_DIR \
        -DCMAKE_INSTALL_PREFIX=$prefix \
        -DPETSC_DIR=$PETSC_DIR \
        -DPETSC_ARCH=$PETSC_ARCH \
        -DAMGX_DIR=$AMGX_DIR \
        -DCUDA_DIR=$CUDA_DIR \
        -DBUILD_SHARED_LIBS=$VAL > /dev/null 2>&1
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

AC_CHECK_HEADER([AmgXSolver.hpp],
                [],
                [AC_MSG_ERROR([Couldn't find AmgXSolver.hpp...
Please use '--with-amgxwrapper-dir=PATH' to provide the directory of the package
or '--enable-amgxwrapper' to download and install AmgXWrapper-1.1 (AmgX-2.0).
])])

AMGXWRAPPER_LIBS="-lAmgXWrapper -lamgxsh -lcudart -lcublas -lcusparse"
AC_SUBST(AMGXWRAPPER_LIBS, $AMGXWRAPPER_LIBS)

PACKAGE_RESTORE_ENVIRONMENT

echo

])
