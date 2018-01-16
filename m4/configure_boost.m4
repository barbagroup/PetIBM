# CONFIGURE_BOOST
# ---------------
# brief: Configures required package boost.

AC_DEFUN([CONFIGURE_BOOST],[

echo
echo "========================="
echo "Configuring Boost library"
echo "========================="

PACKAGE_INITIALIZE_ENVIRONMENT
PACKAGE_SETUP_ENVIRONMENT

AC_ARG_VAR([BOOST_DIR], [Boost directory])

AC_ARG_WITH([boost-dir],
            AS_HELP_STRING([--with-boost-dir=PATH],
                           [set Boost include directory]),
            [BOOST_DIR=$withval],
            [])

AC_ARG_ENABLE([boost],
              AS_HELP_STRING([--enable-boost],
                             [download and install Boost-1.65.1]),
              [IS_INSTALL_BOOST=yes],
              [IS_INSTALL_BOOST=no])

if test "x$BOOST_DIR" != "x"; then
  if test -d "$BOOST_DIR"; then
    SETUP_BOOST
  else
    AC_MSG_ERROR([
$BOOST_DIR is not a valid directory;
please use '--with-boost-dir=PATH' to provide the directory of the package])
  fi
else
  if test "x$IS_INSTALL_BOOST" = "xyes"; then
    INSTALL_BOOST
  fi
fi

AC_CHECK_HEADER([boost/shared_ptr.hpp],
                [],
                [AC_MSG_ERROR([Couldn't find boost/shared_ptr.hpp...
Please use '--with-boost-dir=PATH' to provide the directory of the package
or '--enable-boost' to download and use Boost-1.65.1.
])])

PACKAGE_RESTORE_ENVIRONMENT

echo

]) # CONFIGURE_BOOST


AC_DEFUN([SETUP_BOOST], [

AC_SUBST(BOOST_DIR, $BOOST_DIR)
AC_MSG_NOTICE([using Boost: $BOOST_DIR])
BOOST_CPPFLAGS="-I$BOOST_DIR"
AC_SUBST(BOOST_CPPFLAGS, $BOOST_CPPFLAGS)
CPPFLAGS_APPEND($BOOST_CPPFLAGS)

]) # SETUP_BOOST


AC_DEFUN([INSTALL_BOOST],[

PACKAGE_INITIALIZE_ENVIRONMENT
PACKAGE_SETUP_ENVIRONMENT
AC_MSG_NOTICE([

 ****************************************
 *             Downloading              *
 *            Boost-1.65.1              *
 ****************************************
])
echo "*** INFO *** Downloading Boost-1.65.1... "
VERSION=1.65.1
TARBALL=boost_1_65_1.tar.gz
URL=https://dl.bintray.com/boostorg/release/1.65.1/source/$TARBALL
wget $URL -P /tmp
BOOST_DIR=$BUILDDIR/externalpackages/boost-$VERSION
AC_SUBST([BOOST_DIR], [$BOOST_DIR])
mkdir -p $BOOST_DIR
tar -xzf /tmp/$TARBALL -C $BOOST_DIR --strip-components=1
rm -f /tmp/$TARBALL
echo "*** INFO *** Installing Boost-1.65.1 headers... "
mkdir -p $prefix/include
cp -r $BOOST_DIR/boost $prefix/include/.
PACKAGE_RESTORE_ENVIRONMENT

echo

]) # INSTALL_BOOST
