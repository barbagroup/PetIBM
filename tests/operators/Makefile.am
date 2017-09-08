check_PROGRAMS = createbnhead-test

AM_CPPFLAGS = -I$(top_srcdir)/include \
              -I$(top_srcdir)/external/gtest-1.8.0/include \
              -I$(top_srcdir)/external/yaml-cpp-0.5.1/include \
              -I$(top_srcdir)/external/boost-1.57.0

LDADD = $(top_builddir)/src/libpetibm.la \
        $(top_builddir)/external/gtest-1.8.0/libgtest.a \
        -lstdc++fs

createbnhead_test_SOURCES = createbnhead_test.cpp
createbnhead_test_CPPFLAGS = $(AM_CPPFLAGS)
createbnhead_test_LDADD = $(LDADD)
