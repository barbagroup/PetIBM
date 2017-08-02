TESTS = delta-test \
        singlebody-test \
        cartesianmesh-test

check_PROGRAMS = delta-test \
                 singlebody-test \
                 cartesianmesh-test

AM_CPPFLAGS = -I$(top_srcdir)/src \
              -I$(top_srcdir)/external/gtest-1.8.0/include \
              -I$(top_srcdir)/external/yaml-cpp-0.5.1/include

LDADD = $(top_builddir)/src/utilities/libutilities.a \
        $(top_builddir)/external/yaml-cpp-0.5.1/libyamlcpp.a \
        $(top_builddir)/external/gtest-1.8.0/libgtest.a \
        -lstdc++fs

delta_test_DEPENDENCIES = input_data
delta_test_SOURCES = delta_test.cpp
delta_test_CPPFLAGS = $(AM_CPPFLAGS)
delta_test_LDADD = $(LDADD)

singlebody_test_DEPENDENCIES = input_data
singlebody_test_SOURCES = singlebody_test.cpp
singlebody_test_CPPFLAGS = $(AM_CPPFLAGS)
singlebody_test_LDADD = $(LDADD)

cartesianmesh_test_DEPENDENCIES = input_data
cartesianmesh_test_SOURCES = cartesianmesh_test.cpp
cartesianmesh_test_CPPFLAGS = $(AM_CPPFLAGS)
cartesianmesh_test_LDADD = $(LDADD)

input_data:
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
		cp -r $(top_srcdir)/tests/utilities/data $(PWD) ; \
	fi ;

.PHONY: input_data
