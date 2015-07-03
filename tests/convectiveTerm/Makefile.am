bin_PROGRAMS = convectiveTermTest2d \
               convectiveTermTest3d

convectiveTermTest2d_SOURCES = ConvectiveTerm.cpp \
                               convectiveTermTest.cpp

convectiveTermTest2d_CPPFLAGS = -I$(top_srcdir)/src/utilities \
                                -I$(top_srcdir)/src/solvers \
                                -D DIMENSIONS=2

convectiveTermTest2d_CXXFLAGS = -std=c++0x

convectiveTermTest2d_LDADD = $(top_srcdir)/src/utilities/libutilities.a \
                             $(top_srcdir)/src/solvers/libsolvers.a

convectiveTermTest3d_SOURCES = ConvectiveTerm.cpp \
                               convectiveTermTest.cpp

convectiveTermTest3d_CPPFLAGS = -I$(top_srcdir)/src/utilities \
                                -I$(top_srcdir)/src/solvers \
                                -D DIMENSIONS=3

convectiveTermTest3d_CXXFLAGS = -std=c++0x

convectiveTermTest3d_LDADD = $(top_srcdir)/src/utilities/libutilities.a \
                             $(top_srcdir)/src/solvers/libsolvers.a