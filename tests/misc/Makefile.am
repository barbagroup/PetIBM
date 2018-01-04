check_PROGRAMS = delta-test

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

delta_test_SOURCES = delta_test.cpp
delta_test_CPPFLAGS = $(AM_CPPFLAGS)
delta_test_LDADD = $(LADD)
