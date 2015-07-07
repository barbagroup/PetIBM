# CONFIGURE_BOOST
# ---------------
# brief: Configures required package Boost.
AC_DEFUN([CONFIGURE_BOOST], [

echo
echo "=================================="
echo "Configuring required package Boost"
echo "=================================="

PACKAGE_SETUP_ENVIRONMENT

AC_ARG_VAR(BOOST_ROOT, [location of Boost installation.])

echo "check for a system Boost library"
$as_unset boost_cv_inc_path
BOOST_REQUIRE([1.55], [AC_MSG_ERROR([could not find system Boost library])])

BOOST_ARRAY
BOOST_MULTIARRAY

PACKAGE_CPPFLAGS_PREPEND("$BOOST_CPPFLAGS")

PACKAGE_RESTORE_ENVIRONMENT

])