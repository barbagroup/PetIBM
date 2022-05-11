# CONFIGURE_YAMLCPP
# -----------------
# brief: Configures required package yaml-cpp.
AC_DEFUN([CONFIGURE_YAMLCPP], [

echo
echo "====================================="
echo "Configuring required package yaml-cpp"
echo "====================================="

PACKAGE_INITIALIZE_ENVIRONMENT
PACKAGE_SETUP_ENVIRONMENT

AC_ARG_WITH([yamlcpp-dir],
            AS_HELP_STRING([--with-yamlcpp-dir=PATH],
                           [set yaml-cpp directory]),
            [YAMLCPP_DIR=$withval],
            [])

AC_ARG_WITH([yamlcpp-include],
            AS_HELP_STRING([--with-yamlcpp-include=PATH],
                           [set yaml-cpp include directory]),
            [YAMLCPP_INCLUDE_DIR=$withval],
            [])

AC_ARG_WITH([yamlcpp-lib],
            AS_HELP_STRING([--with-yamlcpp-lib=PATH],
                           [set yaml-cpp lib directory]),
            [YAMLCPP_LIB_DIR=$withval],
            [])

AC_ARG_ENABLE([yamlcpp],
              AS_HELP_STRING([--enable-yamlcpp],
                             [download and install yaml-cpp-0.6.2]),
              [IS_INSTALL_YAMLCPP=yes],
              [IS_INSTALL_YAMLCPP=no])

if test "x$YAMLCPP_DIR" != "x"; then
  CHECK_YAMLCPP_DIR($YAMLCPP_DIR)
else
  if test "x$YAMLCPP_INCLUDE_DIR" != "x"; then
    CHECK_YAMLCPP_INCLUDE_DIR($YAMLCPP_INCLUDE_DIR)
  fi
  if test "x$YAMLCPP_LIB_DIR" != "x"; then
    CHECK_YAMLCPP_LIB_DIR($YAMLCPP_LIB_DIR)
  fi
  if test "x$IS_INSTALL_YAMLCPP" = "xyes"; then
    INSTALL_YAMLCPP
  fi
fi

AC_CHECK_HEADER([yaml-cpp/yaml.h], [],
                [AC_MSG_ERROR([could not find yaml-cpp/yaml.h...
Use '--with-yamlcpp-dir=PATH' to provide the directory of the package or
Use '--with-yamlcpp-include=PATH' and '--with-yamlcpp-lib=PATH' or
Use '--enable-yamlcpp' to download, build, and install yaml-cpp-0.6.2.])])

AC_SUBST(YAMLCPP_LIBS, -lyaml-cpp)

PACKAGE_RESTORE_ENVIRONMENT

echo

]) # CONFIGURE_YAMLCPP


# CHECK_YAMLCPP_DIR
# -----------------
# brief: Checks the existence of the yaml-cpp directory.
AC_DEFUN([CHECK_YAMLCPP_DIR], [

YAMLCPP_DIR=$1
AC_MSG_NOTICE([using yaml-cpp: $YAMLCPP_DIR])

if test ! -d "$YAMLCPP_DIR"; then
  AC_MSG_ERROR([$YAMLCPP_DIR is not a valid directory...
Use '--with-yamlcpp-dir=PATH' to provide the directory of the package.])
fi
CHECK_YAMLCPP_INCLUDE_DIR($YAMLCPP_DIR/include)
if test -d "$YAMLCPP_DIR/lib64"; then
  CHECK_YAMLCPP_LIB_DIR($YAMLCPP_DIR/lib64)
else
  CHECK_YAMLCPP_LIB_DIR($YAMLCPP_DIR/lib)
fi

]) # CHECK_YAMLCPP_DIR


# CHECK_YAMLCPP_INCLUDE_DIR
# -------------------------
# brief: Checks the existence of the yaml-cpp include directory.
AC_DEFUN([CHECK_YAMLCPP_INCLUDE_DIR], [

YAMLCPP_INCLUDE_DIR=$1
AC_MSG_NOTICE([using yaml-cpp include directory: $YAMLCPP_INCLUDE_DIR])

if test -d "$YAMLCPP_INCLUDE_DIR"; then
      YAMLCPP_CPPFLAGS="-I$YAMLCPP_INCLUDE_DIR"
      AC_SUBST(YAMLCPP_CPPFLAGS, $YAMLCPP_CPPFLAGS)
      CPPFLAGS_APPEND($YAMLCPP_CPPFLAGS)
else
  AC_MSG_ERROR([$YAMLCPP_INCLUDE_DIR is not a valid directory...
please use '--with-yamlcpp-include=PATH' to provide the include directory of the package.])
fi

]) # CHECK_YAMLCPP_INCLUDE_DIR


# CHECK_YAMLCPP_LIB_DIR
# ---------------------
# brief: Checks the existence of the yaml-cpp lib directory.
AC_DEFUN([CHECK_YAMLCPP_LIB_DIR], [

YAMLCPP_LIB_DIR=$1
AC_MSG_NOTICE([using yaml-cpp lib directory: $YAMLCPP_LIB_DIR])

if test -d "$YAMLCPP_LIB_DIR"; then
      YAMLCPP_LDFLAGS="-L$YAMLCPP_LIB_DIR -Wl,-rpath,$YAMLCPP_LIB_DIR"
      AC_SUBST(YAMLCPP_LDFLAGS, $YAMLCPP_LDFLAGS)
      LDFLAGS_APPEND($YAMLCPP_LDFLAGS)
else
  AC_MSG_ERROR([$YAMLCPP_LIB_DIR is not a valid directory...
Use '--with-yamlcpp-lib=PATH' to provide the lib directory of the package.])
fi

]) # CHECK_YAMLCPP_LIB_DIR


# INSTALL_YAMLCPP
# ---------------
# brief: Downloads, builds, and installs yaml-cpp.
AC_DEFUN([INSTALL_YAMLCPP], [

version=0.6.2
tarball=yaml-cpp-$version.tar.gz
url=https://github.com/jbeder/yaml-cpp/archive/$tarball
YAMLCPP_DIR=$BUILDDIR/externalpackages/yaml-cpp-$version
if test ! -d "$YAMLCPP_DIR"; then
  echo "downloading yaml-cpp-$version... "
  wget -q $url -P /tmp
  mkdir -p $YAMLCPP_DIR/build
  tar -xzf /tmp/$tarball -C $YAMLCPP_DIR --strip-components=1
  rm -f /tmp/$tarball
fi
echo "building and installing yaml-cpp-$version... "
cd $YAMLCPP_DIR/build
if test "x$enable_shared" = "xyes"; then
  cmake $YAMLCPP_DIR \
    -DCMAKE_BUILD_TYPE="Release" \
    -DCMAKE_INSTALL_PREFIX=$prefix \
    -DBUILD_SHARED_LIBS=ON \
    -DYAML_CPP_BUILD_TESTS=OFF \
    -DCMAKE_MACOSX_RPATH=1
  make all -j
  make install
fi
if test "x$enable_static" = "xyes"; then
  cmake $YAMLCPP_DIR \
    -DCMAKE_BUILD_TYPE="Release" \
    -DCMAKE_INSTALL_PREFIX=$prefix \
    -DBUILD_SHARED_LIBS=OFF \
    -DYAML_CPP_BUILD_TESTS=OFF \
    -DCMAKE_MACOSX_RPATH=1
  make all -j
  make install
fi
cd $BUILDDIR

echo

]) # INSTALL_YAMLCPP
