# CONFIGURE_AMGX
# ---------------
# brief: Configures required package AmgX.

AC_DEFUN([CONFIGURE_AMGX],[

echo
echo "================================="
echo "Configuring optional package AmgX"
echo "================================="

PACKAGE_SETUP_ENVIRONMENT

# check for presence of `--with-amgx=PATH`
AC_ARG_WITH([amgx],
            AS_HELP_STRING([--with-amgx=PATH],
                           [set AmgX directory]),
            [if test ! -d "$withval" ; then
               AC_MSG_ERROR([it is necessary to specify an existing directory 
                             when using --with-amgx=PATH])
             fi
             WITH_AMGX=yes
             AMGX_DIR=$withval],
            [WITH_AMGX=no])

AC_MSG_NOTICE([using AMGX_DIR: $AMGX_DIR])

AMGX_INC_PATH=-I`find $AMGX_DIR -name amgx_c.h | sed -e 's/\/amgx_c.h//g'`
AMGX_LIB_PATH=-L`find $AMGX_DIR -name libamgxsh.so | sed -e 's/\/libamgxsh.so//g'`
AMGX_LIBRARY=-lamgxsh

CPPFLAGS_PREPEND($AMGX_INC_PATH)
LDFLAGS_PREPEND($AMGX_LIB_PATH)
LIBS_PREPEND($AMGX_LIBRARY)

# check for presence of header file amgx_c.h
AC_CHECK_HEADER([amgx_c.h], ,
                AC_MSG_ERROR([could not find header file amgx_c.h]))
# check for presence of library amgxh
AC_CHECK_LIB([amgxsh],
             [AMGX_initialize], ,
             AC_MSG_ERROR([could not find library amgxsh]))

PACKAGE_CPPFLAGS_PREPEND($AMGX_INC_PATH)
PACKAGE_LDFLAGS_PREPEND($AMGX_LIB_PATH)
PACKAGE_LIBS_PREPEND($AMGX_LIBRARY)

PACKAGE_RESTORE_ENVIRONMENT

AC_SUBST(AMGX_DIR, $AMGX_DIR)

echo

])
