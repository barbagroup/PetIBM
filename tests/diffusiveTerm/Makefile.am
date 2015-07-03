bin_PROGRAMS = diffusiveTermTest2d \
               diffusiveTermTest3d

diffusiveTermTest2d_SOURCES = DiffusiveTerm.cpp \
                              diffusiveTermTest.cpp

diffusiveTermTest2d_CPPFLAGS = -I$(top_srcdir)/src/utilities \
                               -I$(top_srcdir)/src/solvers \
                               -D DIMENSIONS=2

diffusiveTermTest2d_CXXFLAGS = -std=c++0x

diffusiveTermTest2d_LDADD = $(top_srcdir)/src/utilities/libutilities.a \
                            $(top_srcdir)/src/solvers/libsolvers.a

diffusiveTermTest3d_SOURCES = DiffusiveTerm.cpp \
                              diffusiveTermTest.cpp

diffusiveTermTest3d_CPPFLAGS = -I$(top_srcdir)/src/utilities \
                               -I$(top_srcdir)/src/solvers \
                               -D DIMENSIONS=3

diffusiveTermTest3d_CXXFLAGS = -std=c++0x

diffusiveTermTest3d_LDADD = $(top_srcdir)/src/utilities/libutilities.a \
                            $(top_srcdir)/src/solvers/libsolvers.a