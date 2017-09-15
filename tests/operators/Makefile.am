check_PROGRAMS = createbnhead-test

AM_CPPFLAGS = -I$(top_srcdir)/include \
              $(YAMLCPP_CPPFLAGS) \
              $(GTEST_CPPFLAGS)

LADD = $(top_builddir)/src/operators/liboperators.la \
       $(YAMLCPP_LDFLAGS) $(YAMLCPP_LIBS) \
       $(GTEST_LDFLAGS) $(GTEST_LIBS)

createbnhead_test_SOURCES = createbnhead_test.cpp
createbnhead_test_CPPFLAGS = $(AM_CPPFLAGS)
createbnhead_test_LDADD = $(LADD)
