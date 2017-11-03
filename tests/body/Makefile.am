check_PROGRAMS = singlebody-test

AM_CPPFLAGS = -I$(top_srcdir)/include \
              $(PETSC_CPPFLAGS) \
              $(YAMLCPP_CPPFLAGS) \
              $(GTEST_CPPFLAGS)

LADD = $(top_builddir)/src/libpetibm.la \
       $(PETSC_LDFLAGS) $(PETSC_LIBS) \
       $(YAMLCPP_LDFLAGS) $(YAMLCPP_LIBS) \
       $(GTEST_LDFLAGS) $(GTEST_LIBS)

singlebody_test_DEPENDENCIES = input_data
singlebody_test_SOURCES = singlebody_test.cpp
singlebody_test_CPPFLAGS = $(AM_CPPFLAGS)
singlebody_test_LDADD = $(LADD)

input_data:
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
		cp $(top_srcdir)/tests/body/body2d.txt $(PWD) ; \
		cp $(top_srcdir)/tests/body/body3d.txt $(PWD) ; \
	fi ;

.PHONY: input_data
