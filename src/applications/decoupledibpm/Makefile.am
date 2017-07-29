noinst_LIBRARIES = libdecoupledibpm.a

libdecoupledibpm_a_SOURCES = decoupledibpm.cpp

libdecoupledibpm_a_CPPFLAGS = -I$(top_srcdir)/src/utilities \
                              -I$(top_srcdir)/src/operators \
                              -I$(top_srcdir)/src/linSolvers \
                              -I$(top_srcdir)/external/yaml-cpp-0.5.1/include \
                              -I$(top_srcdir)/external/boost-1.57.0


bin_PROGRAMS = petibm-decoupledibpm

petibm_decoupledibpm_SOURCES = main.cpp

petibm_decoupledibpm_CPPFLAGS = -I$(top_srcdir)/src/applications \
																-I$(top_srcdir)/src/linSolvers \
																-I$(top_srcdir)/src/operators \
															  -I$(top_srcdir)/src/utilities \
															  -I$(top_srcdir)/external/yaml-cpp-0.5.1/include \
                                -I$(top_srcdir)/external/boost-1.57.0

petibm_decoupledibpm_LDADD = $(top_builddir)/src/applications/decoupledibpm/libdecoupledibpm.a \
														 $(top_builddir)/src/linSolvers/liblinsolvers.a \
														 $(top_builddir)/src/operators/liboperators.a \
														 $(top_builddir)/src/utilities/libutilities.a \
						                 $(top_builddir)/external/yaml-cpp-0.5.1/libyamlcpp.a \
						                -lstdc++fs

if WITH_AMGXWRAPPER
petibm_decoupledibpm_LDADD += $(top_builddir)/external/AmgXWrapper-1.0-beta2/libamgxwrapper.a
endif
