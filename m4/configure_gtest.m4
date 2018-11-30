# CONFIGURE_YAMLCPP
# -----------------
# brief: Configures required package gtest.
AC_DEFUN([CONFIGURE_GTEST], [

echo
echo "=================================="
echo "Configuring required package gtest"
echo "=================================="

PACKAGE_INITIALIZE_ENVIRONMENT
PACKAGE_SETUP_ENVIRONMENT

AC_ARG_WITH([gtest-dir],
            AS_HELP_STRING([--with-gtest-dir=PATH],
                           [set gtest directory]),
            [GTEST_DIR=$withval],
            [])

AC_ARG_WITH([gtest-include],
            AS_HELP_STRING([--with-gtest-include=PATH],
                           [set gtest include directory]),
            [GTEST_INCLUDE_DIR=$withval],
            [])

AC_ARG_WITH([gtest-lib],
            AS_HELP_STRING([--with-gtest-lib=PATH],
                           [set gtest lib directory]),
            [GTEST_LIB_DIR=$withval],
            [])

AC_ARG_ENABLE([gtest],
              AS_HELP_STRING([--enable-gtest],
                             [download and install gtest-1.7.0]),
              [IS_INSTALL_GTEST=yes],
              [IS_INSTALL_GTEST=no])

if test "x$GTEST_DIR" != "x"; then
  CHECK_GTEST_DIR($GTEST_DIR)
else
  if test "x$GTEST_INCLUDE_DIR" != "x"; then
    CHECK_GTEST_INCLUDE_DIR($GTEST_INCLUDE_DIR)
  fi
  if test "x$GTEST_LIB_DIR" != "x"; then
    CHECK_GTEST_LIB_DIR($GTEST_LIB_DIR)
  fi
  if test "x$IS_INSTALL_GTEST" = "xyes"; then
    INSTALL_GTEST
  fi
fi

AC_CHECK_HEADER([gtest/gtest.h], [],
                [AC_MSG_ERROR([Couldn't find gtest/gtest.h...
Use '--with-gtest-dir=PATH' to provide the directory of the package or
Use '--with-gtest-include=PATH' and '--with-gtest-lib=PATH' or
Use '--enable-gtest' to download, build, and install gtest-1.7.0.])])

AC_SUBST(GTEST_LIBS, -lgtest)

PACKAGE_RESTORE_ENVIRONMENT

echo

]) # CONFIGURE_GTEST


# CHECK_GTEST_DIR
# ---------------
# brief: Checks the existence of the gtest directory.
AC_DEFUN([CHECK_GTEST_DIR], [

GTEST_DIR=$1
AC_MSG_NOTICE([using gtest: $GTEST_DIR])

if test ! -d "$GTEST_DIR"; then
  AC_MSG_ERROR([$GTEST_DIR is not a valid directory...
Use '--with-gtest-dir=PATH' to provide the directory of the package.])
fi
CHECK_GTEST_INCLUDE_DIR($GTEST_DIR/include)
if test -d "$GTEST_DIR/lib64"; then
  CHECK_GTEST_LIB_DIR($GTEST_DIR/lib64)
else
  CHECK_GTEST_LIB_DIR($GTEST_DIR/lib)
fi

]) # CHECK_GTEST_DIR


# CHECK_GTEST_INCLUDE_DIR
# -----------------------
# brief: Checks the existence of the gtest include directory.
AC_DEFUN([CHECK_GTEST_INCLUDE_DIR], [

GTEST_INCLUDE_DIR=$1
AC_MSG_NOTICE([using gtest include directory: $GTEST_INCLUDE_DIR])

if test -d "$GTEST_INCLUDE_DIR"; then
      GTEST_CPPFLAGS="-I$GTEST_INCLUDE_DIR"
      AC_SUBST(GTEST_CPPFLAGS, $GTEST_CPPFLAGS)
      CPPFLAGS_APPEND($GTEST_CPPFLAGS)
else
  AC_MSG_ERROR([$GTEST_INCLUDE_DIR is not a valid directory...
Use '--with-gtest-include=PATH' to provide the include directory of the package.])
fi

]) # CHECK_GTEST_INCLUDE_DIR


# CHECK_GTEST_LIB_DIR
# -------------------
# brief: Checks the existence of the gtest lib directory.
AC_DEFUN([CHECK_GTEST_LIB_DIR], [

GTEST_LIB_DIR=$1
AC_MSG_NOTICE([using gtest lib directory: $GTEST_LIB_DIR])

if test -d "$GTEST_LIB_DIR"; then
      GTEST_LDFLAGS="-L$GTEST_LIB_DIR -Wl,-rpath,$GTEST_LIB_DIR"
      AC_SUBST(GTEST_LDFLAGS, $GTEST_LDFLAGS)
      LDFLAGS_APPEND($GTEST_LDFLAGS)
else
  AC_MSG_ERROR([$GTEST_LIB_DIR is not a valid directory...
Use '--with-gtest-lib=PATH' to provide the lib directory of the package.])
fi

]) # CHECK_GTEST_LIB_DIR


# INSTALL_GTEST
# -------------
# brief: Downloads, builds, and installs gtest.
AC_DEFUN([INSTALL_GTEST], [

version=1.7.0
tarball=release-$version.tar.gz
url=https://github.com/google/googletest/archive/$tarball
GTEST_DIR=$BUILDDIR/externalpackages/gtest-$version
if test ! -d "$GTEST_DIR"; then
  echo "downloading gtest-$version... "
  wget -q $url -P /tmp
  mkdir -p $GTEST_DIR/build
  tar -xzf /tmp/$tarball -C $GTEST_DIR --strip-components=1
  rm -f /tmp/$tarball
fi
echo "building and installing gtest-$version... "
cd $GTEST_DIR/build
if test "x$enable_shared" = "xyes"; then
  cmake $GTEST_DIR \
    -DCMAKE_BUILD_TYPE="Release" \
    -DCMAKE_INSTALL_PREFIX=$prefix \
    -DBUILD_SHARED_LIBS=ON \
    -DCMAKE_MACOSX_RPATH=1
  make all
fi
if test "x$enable_static" = "xyes"; then
  cmake $GTEST_DIR \
    -DCMAKE_BUILD_TYPE="Release" \
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

echo

]) # INSTALL_GTEST
