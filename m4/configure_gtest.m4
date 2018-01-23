# CONFIGURE_YAMLCPP
# ---------------
# brief: Configures required package gtest.

AC_DEFUN([CONFIGURE_GTEST],[

echo
echo "=================================="
echo "Configuring required package gtest"
echo "=================================="

PACKAGE_INITIALIZE_ENVIRONMENT
PACKAGE_SETUP_ENVIRONMENT

AC_ARG_VAR([GTEST_DIR], [gtest directory])

AC_ARG_WITH([gtest-dir],
            AS_HELP_STRING([--with-gtest-dir=PATH],
                           [set gtest directory]),
            [GTEST_DIR=$withval],
            [])

AC_ARG_ENABLE([gtest],
              AS_HELP_STRING([--enable-gtest],
                             [download and install gtest-1.7.0]),
              [IS_INSTALL_GTEST=yes],
              [IS_INSTALL_GTEST=no])

if test "x$GTEST_DIR" != "x"; then
  if test -d "$GTEST_DIR"; then
    SETUP_GTEST
  else
    AC_MSG_ERROR([
$GTEST_DIR is not a valid directory;
please use '--with-gtest-dir=PATH' to provide the directory of the package])
  fi
else
  if test "x$IS_INSTALL_GTEST" = "xyes"; then
    INSTALL_GTEST
  fi
fi

AC_CHECK_HEADER([gtest/gtest.h],
                [],
                [AC_MSG_ERROR([Couldn't find gtest/gtest.h...
Please use '--with-gtest-dir=PATH' to provide the directory of the package
or '--enable-gtest' to download and install gtest-1.7.0.
])])

AC_SUBST(GTEST_LIBS, -lgtest)

PACKAGE_RESTORE_ENVIRONMENT

echo

]) # CONFIGURE_GTEST


AC_DEFUN([SETUP_GTEST], [

AC_SUBST(GTEST_DIR, $GTEST_DIR)
AC_MSG_NOTICE([using gtest: $GTEST_DIR])
if test -d "$GTEST_DIR/include"; then
  GTEST_CPPFLAGS="-I$GTEST_DIR/include"
fi
AC_SUBST(GTEST_CPPFLAGS, $GTEST_CPPFLAGS)
CPPFLAGS_APPEND($GTEST_CPPFLAGS)
if test -d "$GTEST_DIR/lib64"; then
  GTEST_LDFLAGS="-L$GTEST_DIR/lib64 -Wl,-rpath,$GTEST_DIR/lib64"
elif test -d "$GTEST_DIR/lib"; then
  GTEST_LDFLAGS="-L$GTEST_DIR/lib -Wl,-rpath,$GTEST_DIR/lib"
fi
AC_SUBST(GTEST_LDFLAGS, $GTEST_LDFLAGS)
LDFLAGS_APPEND($GTEST_LDFLAGS)

]) # SETUP_GTEST


AC_DEFUN([INSTALL_GTEST], [

PACKAGE_INITIALIZE_ENVIRONMENT
PACKAGE_SETUP_ENVIRONMENT
AC_MSG_NOTICE([

 ****************************************
 *      Downloading and installing      *
 *             gtest-1.7.0              *
 ****************************************
])
echo "*** INFO *** Downloading gtest-1.7.0... "
VERSION=1.7.0
TARBALL=release-$VERSION.tar.gz
URL=https://github.com/google/googletest/archive/$TARBALL
wget $URL -P /tmp
GTEST_DIR=$BUILDDIR/externalpackages/gtest-$VERSION
AC_SUBST([GTEST_DIR], [$GTEST_DIR])
mkdir -p $GTEST_DIR/build
tar -xzf /tmp/$TARBALL -C $GTEST_DIR --strip-components=1
rm -f /tmp/$TARBALL
echo "*** INFO *** Building gtest-1.7.0... "
cd $GTEST_DIR/build
if test "x$enable_shared" = "xyes"; then
  cmake $GTEST_DIR \
    -DCMAKE_INSTALL_PREFIX=$prefix \
    -DBUILD_SHARED_LIBS=ON \
    -DCMAKE_MACOSX_RPATH=1
  make all
fi
if test "x$enable_static" = "xyes"; then
  cmake $GTEST_DIR \
    -DCMAKE_INSTALL_PREFIX=$prefix \
    -DBUILD_SHARED_LIBS=OFF \
    -DCMAKE_MACOSX_RPATH=1
  make all
fi
mkdir -p $prefix/include
cp -r $GTEST_DIR/include/gtest $prefix/include/.
mkdir -p $prefix/lib
cp -r $GTEST_DIR/build/*.a $prefix/lib/.
cp -r $GTEST_DIR/build/*.so $prefix/lib/.
cp -r $GTEST_DIR/build/*.dylib $prefix/lib/.
cd $BUILDDIR
PACKAGE_RESTORE_ENVIRONMENT

echo

]) # INSTALL_GTEST
