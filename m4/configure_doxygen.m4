# CONFIGURE_DOXYGEN
# ---------------
# brief: Configures Doxygen.

AC_DEFUN([CONFIGURE_DOXYGEN], [
DOXYGEN_PATH=$PATH
AC_ARG_WITH([doxygen],
  AS_HELP_STRING([--with-doxygen=PATH], 
                 [manually set directory where the doxygen executable resides to PATH]),
 [if test ! -d "$withval" ; then
    AC_MSG_ERROR([it is necessary to specify an existing directory when using --with-doxygen=PATH])
  fi
  DOXYGEN_DIR=$withval
  DOXYGEN_PATH=$DOXYGEN_DIR$PATH_SEPARATOR$DOXYGEN_PATH])
AC_PATH_PROG(DOXYGEN, doxygen, [], $DOXYGEN_PATH)
if test x$DOXYGEN != x ; then
  HAVE_DOXYGEN=yes
else
  AC_MSG_WARN([doxygen not found])
  HAVE_DOXYGEN=no
  echo "if doxygen is install, specify its location via --with-doxygen=PATH"
fi
AC_SUBST(HAVE_DOXYGEN, $HAVE_DOXYGEN)
AC_SUBST(DOXYGEN, $DOXYGEN)
AC_SUBST(DOXYGEN_DIR, $DOXYGEN_DIR)

]) # CONFIGURE_DOXYGEN
