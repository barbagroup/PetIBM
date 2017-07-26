SUBDIRS = utilities \
					operators \
					linSolvers \
					applications

bin_PROGRAMS = petibm2d-navierstokes \
							 petibm3d-navierstokes

petibm2d_navierstokes_SOURCES = main.cpp

petibm2d_navierstokes_CPPFLAGS = -I$(top_srcdir)/src/utilities \
				  											 -I$(top_srcdir)/src/operators \
					 											 -I$(top_srcdir)/src/linSolvers \
																 -I$(top_srcdir)/src/applications \
																 -I$(top_srcdir)/external/yaml-cpp-0.5.1/include \
																 -D DIMENSIONS=2

petibm2d_navierstokes_LDADD = $(top_builddir)/src/applications/libnavierstokes.a \
														  $(top_builddir)/src/utilities/libutilities.a \
														  $(top_builddir)/src/operators/liboperators.a \
														  $(top_builddir)/src/linSolvers/liblinsolvers.a \
						                  $(top_builddir)/external/yaml-cpp-0.5.1/libyamlcpp.a \
						                  -lstdc++fs

if WITH_AMGXWRAPPER
petibm2d_navierstokes_LDADD += $(top_builddir)/external/AmgXWrapper-1.0-beta2/libamgxwrapper.a
endif

petibm3d_navierstokes_SOURCES = main.cpp

petibm3d_navierstokes_CPPFLAGS = -I$(top_srcdir)/src/utilities \
				  											 -I$(top_srcdir)/src/operators \
																 -I$(top_srcdir)/src/linSolvers \
																 -I$(top_srcdir)/src/applications \
																 -I$(top_srcdir)/external/yaml-cpp-0.5.1/include \
																 -D DIMENSIONS=3

petibm3d_navierstokes_LDADD = $(top_builddir)/src/applications/libnavierstokes.a \
														  $(top_builddir)/src/utilities/libutilities.a \
														  $(top_builddir)/src/operators/liboperators.a \
														  $(top_builddir)/src/linSolvers/liblinsolvers.a \
							                $(top_builddir)/external/yaml-cpp-0.5.1/libyamlcpp.a \
							                -lstdc++fs

if WITH_AMGXWRAPPER
petibm3d_navierstokes_LDADD += $(top_builddir)/external/AmgXWrapper-1.0-beta2/libamgxwrapper.a
endif
