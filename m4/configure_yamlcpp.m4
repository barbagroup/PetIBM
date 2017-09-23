# CONFIGURE_YAMLCPP
# ---------------
# brief: Configures required package yaml-cpp.

AC_DEFUN([CONFIGURE_YAMLCPP],[

CONFIGURE_BOOST

echo
echo "====================================="
echo "Configuring required package yaml-cpp"
echo "====================================="

PACKAGE_INITIALIZE_ENVIRONMENT
PACKAGE_SETUP_ENVIRONMENT

AC_ARG_VAR([YAMLCPP_DIR], [yaml-cpp directory])

AC_ARG_WITH([yamlcpp-dir],
            AS_HELP_STRING([--with-yamlcpp-dir=PATH],
                           [set yaml-cpp directory]),
            [YAMLCPP_DIR=$withval],
            [])

AC_ARG_ENABLE([yamlcpp],
              AS_HELP_STRING([--enable-yamlcpp],
                             [download and install yaml-cpp-0.5.1]),
              [DOWNLOAD_YAMLCPP=yes],
              [DOWNLOAD_YAMLCPP=no])

if test "x$YAMLCPP_DIR" != "x"; then
  if test -d "$YAMLCPP_DIR"; then
    AC_SUBST(YAMLCPP_DIR, $YAMLCPP_DIR)
    AC_MSG_NOTICE([using yaml-cpp: $YAMLCPP_DIR])
    YAMLCPP_CPPFLAGS="-I$YAMLCPP_DIR/include $BOOST_CPPFLAGS"
    AC_SUBST(YAMLCPP_CPPFLAGS, $YAMLCPP_CPPFLAGS)
    CPPFLAGS_APPEND($YAMLCPP_CPPFLAGS)
    if test -d "$YAMLCPP_DIR/lib"; then
      YAMLCPP_LDFLAGS="-L$YAMLCPP_DIR/lib -Wl,-rpath,$YAMLCPP_DIR/lib"
    elif test -d "$YAMLCPP_DIR/lib64"; then
      YAMLCPP_LDFLAGS="-L$YAMLCPP_DIR/lib64 -Wl,-rpath,$YAMLCPP_DIR/lib64"
    fi
    AC_SUBST(YAMLCPP_LDFLAGS, $YAMLCPP_LDFLAGS)
    LDFLAGS_APPEND($YAMLCPP_LDFLAGS)
  else
    AC_MSG_ERROR([
$YAMLCPP_DIR is not a valid directory;
please use '--with-yamlcpp-dir=PATH' to provide the directory of the package])
  fi
else
  if test "x$DOWNLOAD_YAMLCPP" = "xyes"; then
    AC_MSG_NOTICE([
 
 ****************************************
 *      Downloading and installing      *
 *           yaml-cpp-0.5.1             *
 ****************************************
])
    PACKAGE_INITIALIZE_ENVIRONMENT
    PACKAGE_SETUP_ENVIRONMENT
    echo "downloading yaml-cpp-0.5.1... "
    VERSION=0.5.1
    TARBALL=release-$VERSION.tar.gz
    URL=https://github.com/jbeder/yaml-cpp/archive/$TARBALL
    wget $URL -P /tmp
    YAMLCPP_DIR=$BUILDDIR/externalpackages/yaml-cpp-$VERSION
    AC_SUBST([YAMLCPP_DIR], [$YAMLCPP_DIR])
    mkdir -p $YAMLCPP_DIR/build
    tar -xzf /tmp/$TARBALL -C $YAMLCPP_DIR --strip-components=1
    rm -f /tmp/$TARBALL
    echo "building yaml-cpp-0.5.1... "
    cd $YAMLCPP_DIR/build
    for VAL in "ON" "OFF"
    do
      if test ! "x$BOOST_DIR" = "x"; then
        cmake $YAMLCPP_DIR \
          -DCMAKE_INSTALL_PREFIX=$prefix \
          -DBoost_INCLUDE_DIR=$BOOST_DIR \
          -DBUILD_SHARED_LIBS=$VAL \
          -DCMAKE_MACOSX_RPATH=1
      else
        cmake $YAMLCPP_DIR \
          -DCMAKE_INSTALL_PREFIX=$prefix \
          -DBUILD_SHARED_LIBS=$VAL \
          -DCMAKE_MACOSX_RPATH=1
      fi
      make all -j4
      make install
    done
    echo "done! "
    cd $BUILDDIR
    PACKAGE_RESTORE_ENVIRONMENT
  fi
fi

AC_CHECK_HEADER([yaml-cpp/yaml.h],
                [],
                [AC_MSG_ERROR([Couldn't find yaml-cpp/yaml.h...
Please use '--with-yamlcpp-dir=PATH' to provide the directory of the package
or '--enable-yamlcpp' to download and install yaml-cpp-0.5.1.
])])

AC_SUBST(YAMLCPP_LIBS, -lyaml-cpp)

PACKAGE_RESTORE_ENVIRONMENT

echo

])
