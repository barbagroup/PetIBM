check_PROGRAMS = cartesianmesh-test

AM_CPPFLAGS = -I$(top_srcdir)/include \
              $(PETSC_CPPFLAGS) \
              $(YAMLCPP_CPPFLAGS) \
              $(GTEST_CPPFLAGS)

LADD = $(top_builddir)/src/libpetibm.la \
       $(PETSC_LDFLAGS) $(PETSC_LIBS) \
       $(YAMLCPP_LDFLAGS) $(YAMLCPP_LIBS) \
       $(GTEST_LDFLAGS) $(GTEST_LIBS)

cartesianmesh_test_SOURCES = main.cpp \
                             cartesianmesh2d_dirichlet.cpp \
                             cartesianmesh2d_yperiodic.cpp \
                             cartesianmesh3d_dirichlet.cpp
cartesianmesh_test_CPPFLAGS = $(AM_CPPFLAGS)
cartesianmesh_test_LDADD = $(LADD)
