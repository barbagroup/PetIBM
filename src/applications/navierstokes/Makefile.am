noinst_LIBRARIES = libnavierstokes.a

libnavierstokes_a_SOURCES = navierstokes.cpp

libnavierstokes_a_CPPFLAGS = -I$(top_srcdir)/src/utilities \
                             -I$(top_srcdir)/src/operators \
                             -I$(top_srcdir)/src/linSolvers \
                             -I$(top_srcdir)/external/yaml-cpp-0.5.1/include

bin_PROGRAMS = petibm-navierstokes

petibm_navierstokes_SOURCES = main.cpp

petibm_navierstokes_CPPFLAGS = -I$(top_srcdir)/src/utilities \
				  										 -I$(top_srcdir)/src/operators \
					 										 -I$(top_srcdir)/src/linSolvers \
															 -I$(top_srcdir)/src/applications \
															 -I$(top_srcdir)/external/yaml-cpp-0.5.1/include

petibm_navierstokes_LDADD = $(top_builddir)/src/applications/navierstokes/libnavierstokes.a \
														$(top_builddir)/src/utilities/libutilities.a \
														$(top_builddir)/src/operators/liboperators.a \
														$(top_builddir)/src/linSolvers/liblinsolvers.a \
						                $(top_builddir)/external/yaml-cpp-0.5.1/libyamlcpp.a \
						                -lstdc++fs

if WITH_AMGXWRAPPER
petibm_navierstokes_LDADD += $(top_builddir)/external/AmgXWrapper-1.0-beta2/libamgxwrapper.a
endif
