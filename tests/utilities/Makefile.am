check_PROGRAMS = delta-test \
                 singlebody-test \
                 cartesianmesh-test \
                 singleboundary-test

AM_CPPFLAGS = -I$(top_srcdir)/include \
              $(PETSC_CPPFLAGS) \
              $(YAMLCPP_CPPFLAGS) \
              $(GTEST_CPPFLAGS)

LADD = $(top_builddir)/src/utilities/libutilities.la \
       $(PETSC_LDFLAGS) $(PETSC_LIBS) \
       $(YAMLCPP_LDFLAGS) $(YAMLCPP_LIBS) \
       $(GTEST_LDFLAGS) $(GTEST_LIBS)

delta_test_DEPENDENCIES = input_data
delta_test_SOURCES = delta_test.cpp
delta_test_CPPFLAGS = $(AM_CPPFLAGS)
delta_test_LDADD = $(LADD)

singlebody_test_DEPENDENCIES = input_data
singlebody_test_SOURCES = singlebody_test.cpp
singlebody_test_CPPFLAGS = $(AM_CPPFLAGS)
singlebody_test_LDADD = $(LADD)

cartesianmesh_test_DEPENDENCIES = input_data
cartesianmesh_test_SOURCES = cartesianmesh_test.cpp
cartesianmesh_test_CPPFLAGS = $(AM_CPPFLAGS)
cartesianmesh_test_LDADD = $(LADD)

singleboundary_test_DEPENDENCIES = input_data
singleboundary_test_SOURCES = singleboundary_test.cpp
singleboundary_test_CPPFLAGS = $(AM_CPPFLAGS)
singleboundary_test_LDADD = $(LADD)

input_data:
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
		cp -r $(top_srcdir)/tests/utilities/data $(PWD) ; \
	fi ;

.PHONY: input_data
