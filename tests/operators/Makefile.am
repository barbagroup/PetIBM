check_PROGRAMS = createBnHead-test

AM_CPPFLAGS = -I$(top_srcdir)/src \
              -I$(top_srcdir)/external/gtest-1.8.0/include \
              -I$(top_srcdir)/external/yaml-cpp-0.5.1/include \
              -I$(top_srcdir)/external/boost-1.57.0

LDADD = $(top_builddir)/src/libpetibm.la \
        $(top_builddir)/external/gtest-1.8.0/libgtest.a \
        -lstdc++fs

createBnHead_test_SOURCES = createBnHead_test.cpp
createBnHead_test_CPPFLAGS = $(AM_CPPFLAGS)
createBnHead_test_LDADD = $(LDADD)
