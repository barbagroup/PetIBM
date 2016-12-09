SUBDIRS = utilities solvers

bin_PROGRAMS = petibm2d \
               petibm3d

petibm2d_SOURCES = PetIBM.cpp

petibm2d_CPPFLAGS = -I$(srcdir)/utilities \
                    -I$(srcdir)/solvers \
                    -D DIMENSIONS=2

petibm2d_LDADD = solvers/libsolvers.a \
                 utilities/libutilities.a \
                 $(top_builddir)/external/yaml-cpp-0.5.1/libyamlcpp.a

if WITH_AMGXWRAPPER
petibm2d_LDADD += $(top_builddir)/external/AmgXWrapper-1.0-beta/libamgxwrapper.a
endif

petibm3d_SOURCES = PetIBM.cpp

petibm3d_CPPFLAGS = -I$(srcdir)/utilities \
                    -I$(srcdir)/solvers \
                    -D DIMENSIONS=3

petibm3d_LDADD = solvers/libsolvers.a \
                 utilities/libutilities.a \
                 $(top_builddir)/external/yaml-cpp-0.5.1/libyamlcpp.a

if WITH_AMGXWRAPPER
petibm3d_LDADD += $(top_builddir)/external/AmgXWrapper-1.0-beta/libamgxwrapper.a
endif
