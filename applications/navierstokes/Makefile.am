bin_PROGRAMS = petibm-navierstokes

petibm_navierstokes_SOURCES = main.cpp \
                              navierstokes.cpp

petibm_navierstokes_CPPFLAGS = -I$(top_srcdir)/src \
                               -I$(top_srcdir)/external/yaml-cpp-0.5.1/include \
                               -I$(top_srcdir)/external/boost-1.57.0

#petibm_navierstokes_LDADD = $(top_builddir)/src/linSolvers/liblinsolvers.la \
#                            $(top_builddir)/src/operators/liboperators.la \
#                            $(top_builddir)/src/utilities/libutilities.la \
#                            $(top_builddir)/external/yaml-cpp-0.5.1/libyamlcpp.la \
#                            -lstdc++fs

#if WITH_AMGXWRAPPER
#petibm_navierstokes_LDADD += $(top_builddir)/external/AmgXWrapper-1.0-beta2/libamgxwrapper.la
#endif

petibm_navierstokes_LDADD = $(top_builddir)/src/libpetibm.la \
                            -lstdc++fs
