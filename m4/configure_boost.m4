# CONFIGURE_BOOST
# ---------------
# brief: Configures required package boost.

AC_DEFUN([CONFIGURE_BOOST],[

echo
echo "========================="
echo "Configuring Boost library"
echo "========================="

PACKAGE_SETUP_ENVIRONMENT

AC_ARG_VAR([BOOST_DIR], [Boost directory])

AC_ARG_WITH([boost-dir],
            AS_HELP_STRING([--with-boost-dir=PATH],
                           [set Boost directory]),
            [BOOST_DIR=$withval],
            [])

if test ! "x$BOOST_DIR" = "x"; then
  if test ! -d "$BOOST_DIR"; then
    AC_MSG_ERROR([
$BOOST_DIR is not a valid directory;
please use '--with-boost-dir' to provide the directory])
  else
    CPPFLAGS_PREPEND(-I$BOOST_DIR)
  fi
fi
AC_CHECK_HEADER([boost/shared_ptr.hpp], [DOWNLOAD_BOOST=no], [DOWNLOAD_BOOST=yes])

if test "x$DOWNLOAD_BOOST" = "xyes"; then
  AC_MSG_WARN([downloading Boost-1.65.1])
  VERSION=1.65.1
  TARBALL=boost_1_65_1.tar.gz
  URL=https://dl.bintray.com/boostorg/release/1.65.1/source/$TARBALL
  wget $URL -P /tmp
  BOOST_DIR=externalpackages/boost-$VERSION
  mkdir -p $BOOST_DIR
  tar -xzf /tmp/$TARBALL -C $BOOST_DIR --strip-components=1
  rm -f /tmp/$TARBALL
  mkdir -p $prefix/include
  cp -r $BOOST_DIR/boost $prefix/include/.
  CPPFLAGS_PREPEND(-I$prefix/include)
fi

AS_UNSET([ac_cv_header_boost_shared_ptr_hpp])
AC_CHECK_HEADER([boost/shared_ptr.hpp],
                [],
                [AC_MSG_ERROR([Couldn't find boost/shared_ptr.hpp...])])

AC_SUBST(BOOST_DIR, $BOOST_DIR)

echo

])
