bin_PROGRAMS = navierstokes2d \
							 navierstokes3d

navierstokes2d_SOURCES = navierstokes.cpp

navierstokes2d_CPPFLAGS = -I$(top_srcdir)/src/utilities \
													-I$(top_srcdir)/src/operators \
													-I$(top_srcdir)/src/linSolvers \
													-I$(top_srcdir)/external/yaml-cpp-0.5.1/include \
													-D DIMENSIONS=2

navierstokes2d_LDADD = $(top_builddir)/src/utilities/libutilities.a \
											 $(top_builddir)/src/operators/liboperators.a \
											 $(top_builddir)/src/linSolvers/liblinsolvers.a \
			                 $(top_builddir)/external/yaml-cpp-0.5.1/libyamlcpp.a

if WITH_AMGXWRAPPER
navierstokes2d_LDADD += $(top_builddir)/external/AmgXWrapper-1.0-beta2/libamgxwrapper.a
endif

navierstokes3d_SOURCES = navierstokes.cpp

navierstokes3d_CPPFLAGS = -I$(top_srcdir)/src/utilities \
													-I$(top_srcdir)/src/operators \
													-I$(top_srcdir)/src/linSolvers \
													-I$(top_srcdir)/external/yaml-cpp-0.5.1/include \
													-D DIMENSIONS=3

navierstokes3d_LDADD = $(top_builddir)/src/utilities/libutilities.a \
											 $(top_builddir)/src/operators/liboperators.a \
											 $(top_builddir)/src/linSolvers/liblinsolvers.a \
			                 $(top_builddir)/external/yaml-cpp-0.5.1/libyamlcpp.a

if WITH_AMGXWRAPPER
navierstokes3d_LDADD += $(top_builddir)/external/AmgXWrapper-1.0-beta2/libamgxwrapper.a
endif
