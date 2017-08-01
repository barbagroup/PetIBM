TESTS = utilities-tests

check_PROGRAMS = utilities-tests

utilities_tests_DEPENDENCIES = input_data

utilities_tests_SOURCES = main.cpp

utilities_tests_CPPFLAGS = -I$(top_srcdir)/src \
                           -I$(top_srcdir)/external/gtest-1.8.0/include \
                           -I$(top_srcdir)/external/yaml-cpp-0.5.1/include

utilities_tests_LDADD = $(top_builddir)/src/utilities/libutilities.a \
                        $(top_builddir)/external/yaml-cpp-0.5.1/libyamlcpp.a \
                        $(top_builddir)/external/gtest-1.8.0/libgtest.a \
                        -lstdc++fs

input_data:
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
		cp -r $(top_srcdir)/tests/utilities/data $(PWD) ; \
	fi ;

.PHONY: input_data
