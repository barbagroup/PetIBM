check_PROGRAMS = navierStokesTest

navierStokesTest_SOURCES = NavierStokesTest.cpp

navierStokesTest_DEPENDENCIES = input_data

navierStokesTest_CPPFLAGS = -I$(top_srcdir)/src/utilities \
														-I$(top_srcdir)/src/solvers \
														-I$(top_srcdir)/external/gtest-1.7.0/include \
														-I$(top_srcdir)/external/yaml-cpp-0.5.1/include

if WITH_AMGXWRAPPER
navierStokesTest_CPPFLAGS += -I$(top_srcdir)/external/AmgXWrapper-1.0-beta/src
endif

navierStokesTest_LDADD = $(top_builddir)/src/solvers/libsolvers.a \
												 $(top_builddir)/src/utilities/libutilities.a \
												 $(top_builddir)/external/gtest-1.7.0/libgtest.a \
												 $(top_builddir)/external/yaml-cpp-0.5.1/libyamlcpp.a

if WITH_AMGXWRAPPER
navierStokesTest_LDADD += $(top_builddir)/external/AmgXWrapper-1.0-beta/libamgxwrapper.a
endif

input_data:
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
		cp -r $(top_srcdir)/tests/NavierStokes/case $(PWD) ; \
	fi ;

.PHONY: input_data
