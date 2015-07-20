SUBDIRS = utilities solvers

bin_PROGRAMS = petibm2d \
               petibm3d

petibm2d_SOURCES = PetIBM.cpp

petibm2d_CPPFLAGS = -I$(srcdir)/utilities \
                    -I$(srcdir)/solvers \
                    -D DIMENSIONS=2

petibm2d_CXXFLAGS = -std=c++0x

petibm2d_LDADD = utilities/libutilities.a \
                 solvers/libsolvers.a \
                 $(top_builddir)/external/yaml-cpp-0.5.1/libyamlcpp.a

petibm3d_SOURCES = PetIBM.cpp

petibm3d_CPPFLAGS = -I$(srcdir)/utilities \
                    -I$(srcdir)/solvers \
                    -D DIMENSIONS=3

petibm3d_CXXFLAGS = -std=c++0x

petibm3d_LDADD = utilities/libutilities.a \
                 solvers/libsolvers.a \
                 $(top_builddir)/external/yaml-cpp-0.5.1/libyamlcpp.a
