# CONFIGURE_YAMLCPP
# ---------------
# brief: Configures required package yaml-cpp.

AC_DEFUN([CONFIGURE_YAMLCPP],[

CONFIGURE_BOOST

echo
echo "====================================="
echo "Configuring required package yaml-cpp"
echo "====================================="

PACKAGE_SETUP_ENVIRONMENT

AC_ARG_VAR([YAMLCPP_DIR], [yaml-cpp directory])

AC_ARG_WITH([yamlcpp-dir],
            AS_HELP_STRING([--with-yamlcpp-dir=PATH],
                           [set yaml-cpp directory]),
            [YAMLCPP_DIR=$withval],
            [])

if test ! "x$YAMLCPP_DIR" = "x"; then
  if test ! -d "$YAMLCPP_DIR"; then
    AC_MSG_ERROR([
$YAMLCPP_DIR is not a valid directory;
please use '--with-yamlcpp-dir' to provide the directory])
  else
    CPPFLAGS_PREPEND(-I$YAMLCPP_DIR/include)
  fi
fi
AC_CHECK_FILE([yaml-cpp/yaml.h], [DOWNLOAD_YAMLCPP=no], [DOWNLOAD_YAMLCPP=yes])

if test "x$DOWNLOAD_YAMLCPP" = "xyes"; then
  AC_MSG_WARN([downloading and installing yaml-cpp-0.5.1])
  VERSION=0.5.1
  TARBALL=release-${VERSION}.tar.gz
  URL=https://github.com/jbeder/yaml-cpp/archive/${TARBALL}
  wget ${URL} -P /tmp
  YAMLCPP_DIR=externalpackages/yaml-cpp-${VERSION}
  mkdir -p ${YAMLCPP_DIR}/build
  tar -xzf /tmp/${TARBALL} -C ${YAMLCPP_DIR} --strip-components=1
  rm -f /tmp/${TARBALL}
  cd ${YAMLCPP_DIR}/build
  cmake .. -DCMAKE_INSTALL_PREFIX=$prefix
  make all -j"$(nproc)"
  make install
  cd ../../..
  CPPFLAGS_PREPEND(-I$prefix/include)
fi

AC_CHECK_HEADER([yaml-cpp/yaml.h],
                [],
                [AC_MSG_ERROR([Couldn't find yaml-cpp/yaml.h...])])

AC_SUBST(YAMLCPP_DIR, $YAMLCPP_DIR)

echo

])
