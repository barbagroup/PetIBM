check_PROGRAMS = singleboundary-test

AM_CPPFLAGS = -I$(top_srcdir)/include \
              $(PETSC_CPPFLAGS) \
              $(YAMLCPP_CPPFLAGS) \
              $(GTEST_CPPFLAGS)

LADD = $(top_builddir)/src/libpetibm.la \
       $(PETSC_LDFLAGS) $(PETSC_LIBS) \
       $(YAMLCPP_LDFLAGS) $(YAMLCPP_LIBS) \
       $(GTEST_LDFLAGS) $(GTEST_LIBS)
if WITH_AMGX
LADD += $(AMGXWRAPPER_LDFLAGS) $(AMGXWRAPPER_LIBS)
endif

singleboundary_test_SOURCES = singleboundary_test.cpp
singleboundary_test_CPPFLAGS = $(AM_CPPFLAGS)
singleboundary_test_LDADD = $(LADD)
